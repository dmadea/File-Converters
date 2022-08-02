#!/usr/bin/python3

# Converts the Agilent *.dx HPLC files, reads the absorption (UV prefix) and fluorescence (FLD prefix) data
# if available and saves them into csv files in the same folder. Inspired by the 3rd party Aston library
# https://github.com/bovee/Aston. New files share the structure with the old Agilent .DAD format, but all
# of the actual data are saved as float64 type instead of int16 or int32 as in old format

import sys
import os
import numpy as np
import struct
import xml.etree.ElementTree as ET

from hplc_converter_base import read_utf16, save_mat2csv


def _read_data(file, loc, data_scale_factor):
    """Reads the actual data, file is an open binary file, loc is the location of the data in the file,
    data_scale_factor additionally scales the read data matrix."""

    # various shifts from the initial location for different data
    nrec_shift = 375  # number of records, >i
    sample_name_shift = 955  # utf-16 and the first byte is length of the string
    scale_factor_shift = 3182  # scale factor for data, >d, unit is just after this factor
    data_shift = 4193

    space_len = 22  # length of leading bytes before individual spectra, both for FLD and UV data

    # scaling factors
    scale_time = 1 / 60000  # time
    scale_wl = 1 / 20  # wavelength

    # load number of records (it is also specified in XML tree)
    file.seek(loc + nrec_shift)
    nrec, = struct.unpack('>i', file.read(4))

    # load name of sample, not necessary
    file.seek(loc + sample_name_shift)
    sample_name = read_utf16(file)

    # load scaling factor for data
    file.seek(loc + scale_factor_shift)
    scale_fac, = struct.unpack('>d', file.read(8))
    scale_fac *= data_scale_factor
    # load unit of data
    unit = read_utf16(file)

    # load data itself
    file.seek(loc + data_shift)

    times = np.empty(nrec, dtype=np.float64)
    wavelengths = None
    data_mat = None

    for i in range(nrec):
        # in those 22 bytes before the data, there are 4 same bytes for all spectra:
        # e.g. 46 00 EE 07, last two bytes denotes the size of the block, then there
        # are 4 bytes of time the data corresponds to as <i (le int32), then there are
        # 3 x 2 bytes of little endian uint16 that corresponds to start, end and step of
        # the wavelength range, the remaining bytes are the same for all records and IDK
        # what they are
        leading_bytes = file.read(space_len)
        block_size, = struct.unpack('<H', leading_bytes[2:4])
        times[i],  = struct.unpack('<i', leading_bytes[4:8])  # time of measurement
        if wavelengths is None:
            wl_start, wl_end, wl_step = struct.unpack('<HHH', leading_bytes[8:14])
            wavelengths = np.arange(wl_start, wl_end + wl_step, wl_step) * scale_wl
            data_mat = np.empty((nrec, wavelengths.shape[0]), dtype=np.float64)  # create a data matrix for our data

        # if this is not valid for some type of files, the algorithm needs to be rewritten
        # it assumes the block_size is the same for all records, so far it worked...
        assert (block_size - space_len) / 8 == wavelengths.shape[0]

        # read the block of <d (le float64) values and put them into matrix
        data_mat[i, :] = np.frombuffer(file.read(8 * wavelengths.shape[0]), dtype='<d')

    # apply the scale for data
    data_mat *= scale_fac
    times *= scale_time

    return data_mat, times, wavelengths, sample_name, unit


def process_filepath(fpath):
    _dir, fname = os.path.split(fpath)   # get dir and filename
    fname, _ = os.path.splitext(fname)  # get filename without extension

    with open(fpath, 'rb') as f:
        # read the last 35000 bytes (this should contain all XML tree at the end of DX file)
        block_size = 35000
        f.seek(0, os.SEEK_END)
        f.seek(f.tell() - block_size)

        data = f.read(block_size)

        # find the XML tree
        init = bytes('<ACMD xmlns="urn:schemas-agilent-com:acmd20">', 'utf8')
        end = bytes('</ACMD>', 'utf8')

        start_idx = data.find(init)
        end_idx = data.find(end)

        # get the root of XML tree
        xml_data = data[start_idx:end_idx + len(end)]
        root = ET.fromstring(xml_data)

        ssNs = '{urn:schemas-agilent-com:acmd20}'

        # IDs for location of fluorescence and UV data
        fld_ID = None
        uv_ID = None

        # find the IDs
        for signals in root.iter(ssNs + 'Signal'):
            for item in signals.iter(ssNs + 'Description'):
                if item.text == 'FLD1P,FLD: Spectrum':
                    traceID_el = list(signals.iter(ssNs + "TraceId"))[0]
                    fld_ID = traceID_el.text
                    continue

                if item.text == 'DAD1I,DAD: Spectrum':
                    traceID_el = list(signals.iter(ssNs + "TraceId"))[0]
                    uv_ID = traceID_el.text
                    continue

        # ----------- READ and SAVE FLD DATA -------

        if fld_ID is not None:
            # find the locations of fluorescence and UV data (locations are after the XML tree)
            fld_id_idx = data.find(bytes(fld_ID + '.UVPK', 'utf8'))
            # 4 bytes before the ID is location in data file
            FLD_loc, = struct.unpack('<I', data[fld_id_idx - 4: fld_id_idx])

            data_mat, times, wavelengths, _, unit = _read_data(f, FLD_loc, 1e-6)  # an extra scaling factor

            save_mat2csv(os.path.join(_dir, f'FLD_{fname}.csv'), data_mat, times, wavelengths, unit=unit)

        # ----------- READ and SAVE UV DATA -------

        if uv_ID is not None:
            uv_id_idx = data.find(bytes(uv_ID + '.UVPK', 'utf8'))
            UV_loc, = struct.unpack('<I', data[uv_id_idx - 4: uv_id_idx])

            data_mat, times, wavelengths, _, unit = _read_data(f, UV_loc, 1 / 2000)  # an extra scaling factor

            save_mat2csv(os.path.join(_dir, f'UV_{fname}.csv'), data_mat, times, wavelengths, unit=unit)


if __name__ == '__main__':
    # Show file open dialog if no input files specified on command line
    if len(sys.argv) > 1:
        filepaths = sys.argv[1:]
    else:
        from tkinter import Tk, filedialog

        tk_root = Tk()
        tk_root.withdraw()
        filepaths = list(filedialog.askopenfilenames(parent=None, title='Select dx input files...'))
        tk_root.destroy()
        if len(filepaths) < 1:
            print('No input files selected. Exiting.')
            exit(1)

    for fpath in filepaths:
        process_filepath(fpath)

