#!/usr/bin/python3

# Converts the old Agilent *.UV HPLC files, reads the absorption data and aves them into csv files in the 
# same folder. Inspired by the 3rd party Aston library https://github.com/bovee/Aston. It does not work for
# the new HPLC, files exported as old format (*.D Chemstation) are not read properly for some reason.

import sys
import os
import struct
import numpy as np
from numba import njit
from hplc_converter_base import read_utf16, save_mat2csv

@njit
def get_array(i2arr: np.ndarray, n_wls: int):
    """

    :param i2arr: np.int16 array of raw data
    :param n_wls: number of wavelengths
    :return:
    """
    ret_array = np.empty(n_wls, dtype=np.float64)
    idx = 0
    value = 0
    for i in range(n_wls):
        mv = i2arr[idx]
        if mv == -32768:
            # reinterpret the next two int16 values as int32
            value = i2arr[idx + 1:idx + 3].view(np.int32)[0]
            idx += 3
        else:
            # cumulative add otherwise
            value += mv
            idx += 1
        ret_array[i] = value
    return ret_array


# adapted from https://github.com/bovee/Aston/blob/master/aston/tracefile/agilent_uv.py
# modified with numba function
def _read_data(file):

    space_len = 22  # length of leading bytes before individual spectra
    scale_wl = 1 / 20  # wavelength

    file.seek(0x35A)
    sample_name = read_utf16(file)

    file.seek(0xC15)
    yunit = read_utf16(file)

    file.seek(0x116)
    nrec, = struct.unpack('>i', file.read(4))  # number of records (time points)

    # read data scale factor
    file.seek(0xc0d)
    scale_fac, = struct.unpack('>d', file.read(8))

    times = np.empty(nrec, dtype=np.float64)
    wavelengths = None
    data_mat = None

    file.seek(0x1000)  # data starts here
    for i in range(nrec):
        leading_bytes = file.read(space_len)
        block_size, = struct.unpack('<H', leading_bytes[2:4])
        times[i], = struct.unpack('<i', leading_bytes[4:8])  # time of measurement

        wl_start, wl_end, wl_step = struct.unpack('<HHH', leading_bytes[8:14])
        if wavelengths is None:
            wavelengths = np.arange(wl_start, wl_end + wl_step, wl_step) * scale_wl
            data_mat = np.empty((nrec, wavelengths.shape[0]), dtype=np.float64)  # create a data matrix for our data
        else:
            assert (wl_end - wl_start) // wl_step + 1 == wavelengths.shape[0], "invalid file or different format"

        i2arr = np.frombuffer(file.read(block_size - space_len), dtype='<i2')

        data_mat[i, :] = get_array(i2arr, wavelengths.shape[0])

        # old working code
        # v = 0
        # for j in range(wavelengths.shape[0]):
        #     ov = struct.unpack('<h', file.read(2))[0]
        #     if ov == -32768:
        #         v = struct.unpack('<i', file.read(4))[0]
        #     else:
        #         v += ov
        #     data_mat[i, j] = v

    data_mat *= scale_fac  # / 2000
    times /= 60000

    return data_mat, times, wavelengths, yunit, sample_name


def process_filepath(fpath):

    _dir, fname = os.path.split(fpath)   # get dir and filename
    fname, _ = os.path.splitext(fname)  # get filename without extension

    with open(fpath, 'rb') as f:
        data_mat, elution_times, wavelengths, yunit, name = _read_data(f)

    save_mat2csv(os.path.join(_dir, f'UV_{name}.csv'), data_mat, elution_times, wavelengths, unit=yunit)


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

exit(0)

