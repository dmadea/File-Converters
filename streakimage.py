# taken from https://github.com/nicohofdtz/streakimage

from datetime import datetime, date
from enum import Enum
import re
from types import SimpleNamespace
import numpy as np
import pandas as pd
from collections import OrderedDict, namedtuple
from pathlib import Path

import struct
import os
import configparser
from scipy import interpolate

# from scipy import interpolate


class FileType(Enum):
    """Enum for the file types used by HPD-TA

    can be 8 bit, 16 bit, 32 bit or compressed
    """

    BIT8 = 0
    COMPRESSED = 1  # not used by HPD-TA
    BIT16 = 2
    BIT32 = 3


class CameraModel(Enum):
    """Enum for the camera types."""

    C4742_95 = "uvvis"
    C4742_95_12ER = "ir"


class StreakImage:
    """Parses and holds the data of a recorded image

    args:
        file_path: path to the file to be parsed
        verbose: print verbose comments
        correction: apply camera correction
    """

    def __init__(
        self,
        file_path: str,
        verbose: bool = False,
        correction: bool = False,
        bg=None,
        bg_dict: dict = None,
        title: str = None,
    ):

        self.config = configparser.ConfigParser()
        self.parent_dir = Path(os.path.dirname(__file__))
        self.config.read(self.parent_dir / "streakimage.ini")

        # the hptda fields
        self.file_path = file_path
        self.comment_length: int
        self.width: int
        self.height: int
        self.x_offset: int
        self.y_offset: int
        self.file_type: FileType
        self.data: pd.DataFrame
        self.comment_string: str
        self.parameters: SimpleNamespace

        # the custom fields
        self.verbose = verbose
        self.correction = correction
        self.bg: StreakImage = bg
        self.bg_dict = bg_dict
        self.bg_sub_corr: bool
        self.title = title
        self.transient: pd.Series = None

        self.parse_file(file_path)
        if self.verbose:
            print("Comment length:", self.comment_length)
            print("Width:", self.width, "Height:", self.height)
            print("x-offset:", self.x_offset, "y-offset", self.y_offset)
            print("file-type:", self.file_type.name)
            # print("comment:", self.comment_string)
        if self.bg_dict and not self.bg:
            self.bg = self.bg_dict[self.get_bg_id()]
        if self.bg:
            self.apply_bg_subtraction()
        if self.correction:
            self.apply_camera_correction()
        if self.parameters.Streakcamera.Mode == "Operate":
            self.shift_0_to_max()

    def parse_file(self, file_path: str):
        """Read the given file and parse the content to class fields.

        args:
            file_path: Path to the file to be parsed

        """
        with open(file_path, "rb") as file:
            file_content = file.read()
            self.im_chars = file_content[0:2].decode("utf-8")
            self.comment_length = int.from_bytes(
                file_content[2:4], byteorder="little", signed=False
            )
            self.width = int.from_bytes(
                file_content[4:6], byteorder="little", signed=False
            )
            self.height = int.from_bytes(
                file_content[6:8], byteorder="little", signed=False
            )
            self.x_offset = int.from_bytes(
                file_content[8:9], byteorder="little", signed=False
            )
            self.y_offset = int.from_bytes(
                file_content[10:11], byteorder="little", signed=False
            )
            file_type_int = int.from_bytes(
                file_content[12:13], byteorder="little", signed=False
            )
            self.file_type = FileType(file_type_int)

            nnn = 64 + self.comment_length
            self.comment_string = file_content[64:nnn].decode("utf-8")
            self.parameters = self.parse_comment(self.comment_string)
            self.data = self.parse_data(file_content[nnn:])
            self.bg_sub_corr = self.parameters.Acquisition.BacksubCorr == "1"

    def parse_data(self, binary_data: bytes) -> pd.DataFrame:
        """Parse given data bytes and compile axes

        args:
            binary_data: The part of the file containing the intensity data.

        """
        data = np.zeros((self.height, self.width))
        byte_per_pixel = 2 if self.file_type == FileType.BIT16 else 4
        from_ = 0
        to = byte_per_pixel
        for time in range(0, self.height):
            for wl in range(0, self.width):
                data[time][wl] = int.from_bytes(
                    binary_data[from_:to], byteorder="little", signed=False
                )
                from_ += byte_per_pixel
                to += byte_per_pixel
        wl_axis, time_axis = self.get_axes(binary_data)
        return pd.DataFrame(data=data[::, ::-1], index=time_axis, columns=wl_axis[::-1])

    def get_axes(self, binary_data: bytes) -> tuple:
        Axes = namedtuple("Axes", "wl time")
        if self.parameters.Streakcamera.Mode == "Focus":
            offset = self.width * 4
            wl_axis = np.asarray(
                (struct.unpack_from("f" * (self.width), binary_data[-offset:]))
            )
            time_axis = np.arange(-self.height / 2 + 1, self.height / 2 + 1, 1)
        elif self.parameters.Streakcamera.Mode == "Operate":
            offset = self.width * 4 + self.height * 4
            wl_axis = np.asarray(
                (struct.unpack_from("f" * (self.width), binary_data[-offset:]))
            )
            time_axis = np.asarray(
                (
                    struct.unpack_from(
                        "f" * (self.height), binary_data[(-offset + self.width * 4) :]
                    )
                )
            )
        return Axes(wl=wl_axis, time=time_axis)

    def get_namespace(self, dict_: dict) -> SimpleNamespace:
        return SimpleNamespace(
            **{
                k: (v if not isinstance(v, dict) else self.get_namespace(v))
                for k, v in dict_.items()
            }
        )

    def parse_comment(self, comment: str) -> SimpleNamespace:
        """Parse the comment string and return it as a dictionary

        args:
            comment: the comment string

        return:
            The returned object is a dictionary of dictionaries. The first level is a dictionary of the comment categories.
            The second level dictionaries contain the key-value pairs of each category.

        """
        # build 1st level dictionary
        comment_dict: OrderedDict = OrderedDict()
        # split comment string into category substrings
        category_rex = r"(\[[\w ]+?\],(?:.|\n)*?(?=\[|\Z))"
        comment_list = re.split(category_rex, comment)
        comment_list = [item for item in comment_list if item != ""]

        # build 2nd level dictionaries
        for category in comment_list:
            # remove linebreaks
            category = category.replace("\r", " - ")
            category = category.replace("\n", " - ")
            # category name is written in brackets and is extracted first
            category_match = re.match(r"\[(.*?)\]\,(.*)", category)
            if category_match:
                category_name, category_body = category_match.groups()
            else:
                raise ValueError(
                    f"Category name and/or body could not be parsed.\ncategory string:\n{category}"
                )

            key_rex = r"[a-zA-Z0-9. ]*"
            value_rex = r"[a-zA-Z0-9\- ]+"
            quoted_val_rex = r'"(?:.|\n)*?"'
            comma_or_eos_rex = r"?:,|$| $|"
            key_val_pair_rex = (
                rf"({key_rex})=({value_rex}|{quoted_val_rex})({comma_or_eos_rex})"
            )

            category_dict_with_spaces = OrderedDict(
                re.findall(key_val_pair_rex, category_body)
            )
            category_dict = {
                k.replace(" ", "").replace(".", ""): v.replace('"', "").strip(" -")
                for k, v in category_dict_with_spaces.items()
            }

            comment_dict[category_name.replace(" ", "")] = category_dict

        if self.verbose:
            print("This is the comment string:")
            print(comment)
            print("Comment parsed. This is the resulting dict:")
            with open("params.txt", "w") as params:
                for category in comment_dict:
                    # ToDo: remove
                    # params.write(category + "=namedtuple('" + category + "','")
                    print("-" * 60 + "\n")
                    print(category + ":")
                    print("-" * 60)
                    if category == "comment":
                        print(comment_dict["Comment"]["UserComment"])
                    else:
                        keys: str = ""
                        for key in comment_dict[category]:
                            value = comment_dict[category][key]
                            print("|\t{:.<22s}:{: <28s}".format(key, value) + "|")
                            keys = keys + key + " "
                print("-" * 60 + "\n")
            print(comment_dict)
        param_list = self.get_namespace(comment_dict)
        return param_list

    def is_compatible(self, other: "StreakImage") -> bool:
        """Compare file type and dimension of given classes"""

        if self.height != other.height:
            raise IndexError(
                "Height differs: {:s} vs {:s}".format(self.height, other.height)
            )

        if self.width != other.width:
            raise IndexError(
                "Width differs: {:s} vs {:s}".format(self.width, other.width)
            )

        if self.file_type != other.file_type:
            raise ValueError("File types do not match.")

        return True

    def apply_bg_subtraction(self):
        if self.is_compatible(self.bg):
            self.data = self.data - self.bg.data.values
            self.bg_sub_corr = True
        else:
            raise TypeError("Background is not compatible to data.")

    def apply_manual_offset(self, range_: list):
        a = range_[0][0]
        b = range_[0][1]
        c = range_[1][0]
        d = range_[1][1]

        offset = (self.data.iloc[c:d].loc[:, b:a]).mean().mean()
        self.data -= offset

    def quadratic(self, x, a, pow, c):
        return a * x**pow + c

    def apply_gain_correction(self):
        gain = self.parameters.Streakcamera.MCPGain
        if int(gain) % 10 == 0:
            gcoef = float(self.config["Gain-Correction"][gain])
        else:
            a = float(self.config["Fitted-Gain-Correction"]["a"])
            pow = float(self.config["Fitted-Gain-Correction"]["pow"])
            c = float(self.config["Fitted-Gain-Correction"]["c"])
            gcoef = self.quadratic(float(gain), a, pow, c)
        self.data /= gcoef

    def apply_cornerbg(self, ci1=10, ci2=-10, t_max=10):
        corner_1 = self.data.iloc[:t_max, :ci1]
        corner_2 = self.data.iloc[:t_max, ci2:]
        o1 = corner_1.mean().mean()
        o2 = corner_2.mean().mean()
        offset = (o1 + o2) / 2
        self.data -= offset

    def apply_exp_correction(self):
        suffix_dic = {"ms": 1, "us": 0.001}
        time_and_unit = self.parameters.Acquisition.ExposureTime.split(" ")
        exp_time = int(time_and_unit[0]) * suffix_dic[time_and_unit[1]]
        nr_exp = int(self.parameters.Acquisition.NrExposure)
        cfak = exp_time * nr_exp
        self.data /= cfak

    def get_cam_corr_prefix(self) -> str:
        camera_name = self.parameters.Camera.CameraName
        """Returns the prefix for the cam correction file"""
        date = self.get_date()
        for key in self.config[f"Cam-Correction-Prefix-{camera_name}"]:
            if key != "default":
                start, end = key.split("-")
                if start <= date and date <= end:
                    return self.config[f"Cam-Correction-Prefix-{camera_name}"][key]
        return self.config[f"Cam-Correction-Prefix-{camera_name}"]["default"]

    def apply_camera_correction(self, prefix=None, verbose=False):
        timerange: str = self.parameters.Streakcamera.TimeRange
        prefix = f"_{self.get_cam_corr_prefix()}" if not prefix else f"_{prefix}"
        camera_name = self.parameters.Camera.CameraName
        cam_corrections_folder = self.config["Cam-Corrections-Folders"][camera_name]
        cam_corrections_filename = f"{camera_name}{prefix}_ST{timerange}_correction_{str(self.width)}x{str(self.height)}.npy"
        if verbose:
            print(cam_corrections_filename)

        correction = np.load(
            f"{os.path.dirname(__file__)}/{cam_corrections_folder}/{cam_corrections_filename}"
        )
        self.data = self.data / correction

    def apply_mono_correction(self, extrapolate=False):
        camera_name = self.parameters.Camera.CameraName
        disp = pd.read_csv(
            os.path.join(
                os.path.dirname(__file__),
                f"correction_data/{camera_name}_corrections/{camera_name}_mono.dat",
            ),
            delimiter="\t",
            index_col=0,
            names=[""],
        ).squeeze("columns")
        wls = self.data.columns.values
        if (
            wls[0] < disp.index.values[0] or wls[-1] > disp.index.values[-1]
        ) and not extrapolate:
            raise IndexError(
                "Data bounds exceed the measured mono correction. To use interpolation explicitly set extrapolate=True"
            )
        else:
            interpolate_disp = interpolate.interp1d(
                disp.index.values, disp.values, fill_value="extrapolate"
            )
            corr = interpolate_disp(wls)
            self.data /= corr

    def shift_0_to_max(self):
        self.data.index -= self.index_of_max()

    def index_of_max(self):
        trans = np.sum(self.data, axis=1)
        return trans.idxmax()

    def shift_time_scale(self, shift_value):
        self.data.index += shift_value

    def get_date(self) -> str:
        date_re = re.match(
            r"(?P<month>[0-9]{2})\-(?P<day>[0-9]{2})\-(?P<year>[0-9]{4})",
            self.parameters.Application.Date,
        )
        time_re = re.match(
            r"(?P<hour>[0-9]{2})\:(?P<minute>[0-9]{2}):(?P<second>[0-9]{2})",
            self.parameters.Application.Time,
        )
        if date_re and time_re:
            date = date_re.groupdict()
            time = time_re.groupdict()
        else:
            raise ValueError("Date and/or time could not be parsed.")

        datetime_ = datetime(
            int(date["year"]),
            int(date["month"]),
            int(date["day"]),
            int(time["hour"]),
            int(time["minute"]),
            int(time["second"]),
        )

        return datetime_.date().isoformat().replace("-", "")

    def get_bg_id(self) -> str:
        id = f"bg_ST{self.parameters.Streakcamera.TimeRange}_"
        id += f"g{self.parameters.Streakcamera.MCPGain}_"
        id += f"{self.parameters.Acquisition.NrExposure}x"
        id += f"{self.parameters.Acquisition.ExposureTime.replace(' ','')}"
        return id

    def get_id(self) -> str:
        return os.path.basename(self.file_path)[:12]

    def get_important_params(self, sep=",") -> str:
        spec_keys = f"Grating{sep}Wavelength{sep}SlitWidth"
        spect_values = f"{self.parameters.Spectrograph.Grating}{sep}{self.parameters.Spectrograph.Wavelength}{sep}{self.parameters.Spectrograph.SlitWidth}"

        cam_keys = f"Gain{sep}TimeRange{sep}Mode"
        cam_values = f"{self.parameters.Streakcamera.MCPGain}{sep}{self.parameters.Streakcamera.TimeRange}{sep}{self.parameters.Streakcamera.Mode}"

        acq_keys = f"NrExposure{sep}ExposureTime"
        acq_values = f"{self.parameters.Acquisition.NrExposure}{sep}{self.parameters.Acquisition.ExposureTime}"

        return f"Spectograph\r{spec_keys}\r{spect_values}\rStreakcamera\r{cam_keys}\r{cam_values}\rAcquisition\r{acq_keys}\r{acq_values}"

    def export_wexf(self, path: Path):
        """Export img to the wavelength explicit format readable by Glotaran

        See https://github.com/nicohofdtz/glotaran.github.io/blob/patch-1/legacy/file_formats.md for more information
        """
        with open(
            path,
            "w+",
        ) as file:
            file.write(f"\n\nWavelength explicit\nIntervalnr {self.width}\n")
            self.data.to_csv(file, sep=" ")