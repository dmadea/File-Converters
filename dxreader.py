#!/usr/bin/python3

# Converts the Agilent *.UV files to csv, uses 3rd party Aston library https://github.com/bovee/Aston
# works on files form old HPLC
# does not work for new HPLC, files exported as old format (*.D Chemstation) are not read properly for some reason


import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import struct


# from https://github.com/bovee/Aston
from aston.tracefile import TraceFile

def utf16_read(f):
    """modified string read method which works with our UV files"""
    # determine length to read
    read_len = struct.unpack('>B', f.read(1))[0]
    # read values, decode, and strip
    bytes = f.read(2 * read_len)
    str = bytes.decode('utf16').strip()

    return str


def save_mat2csv(fname, matrix, times=None, wls=None):
    delimiter = ','
    times = np.arange(0, matrix.shape[0]) if times is None else times
    wls = np.arange(0, matrix.shape[1]) if wls is None else wls

    mat = np.hstack((times[:, None], matrix))
    buffer = f'Elution Time | Wavelength->'
    buffer += delimiter + delimiter.join(f"{num}" for num in wls) + '\n'
    buffer += '\n'.join(delimiter.join(f"{num}" for num in row) for row in mat)

    with open(fname, 'w', encoding='utf8') as f:
        f.write(buffer)


def save_bytes(fname, data):
    with open(fname, 'rb+') as f:
        f.write(data)

fname = r'C:\Users\Dominik\Documents\MUNI\Organic Photochemistry\RealTimeSync\Projects\2020-Bilirubin - 2nd half\HPLC\new HPLC\DMbr046, 2Z after workup.sirslt\DMbr046, 2Z after workup.dx'

f = open(fname, 'rb')
f.seek(0x3BB)
str0 = utf16_read(f)
f.seek(0x10D6)
str1 = utf16_read(f)

#  float64 data
# FLD1E,Ambient Temperature
#  <TraceId>25758435-e603-499a-b296-e81c726c10d2</TraceId>

# f.seek(0x1861)
# n_records = int((0x495c1 - 0x1861) / 8)
#
# data = np.frombuffer(f.read(8 * n_records), dtype='<d')  # read little endian float64 values
# data = data.reshape(n_records // 2, 2)
#
# print(data)
# print(data.shape)
#
# diffs = np.diff(data[:, 0])   # save difference across all data

#  FLD1P,FLD: Spectrum
#  <TraceId>083cecfc-187d-4eed-9557-948967ac9c41</TraceId>
#  <ScaleFactor>1E-06</ScaleFactor>
#  <Units>LU</Units>
#  <NumberOfRecords>1306</NumberOfRecords>
space_len = 22
nrec = 1306
scale_fac = 1e-6
scale_time = 1 / 60000
wl_scale = 1 / 20

f.seek(0x4A622)

times = np.empty(nrec, dtype=np.float64)
wls = None
data_mat = None
leading_bytes = None

for i in range(nrec):
    # in those 22 bytes before the data, there are 4 same bytes for all spectra: 46 00 EE 07
    # then there are 4 bytes of time the data corresponds to as <i (le int32)
    # then there are 3 x 2 bytes of little endian uint16 that corresponds to start, end and step
    # of the wavelength range
    # the remaining bytes are the same for all spectra and IDK what they are
    leading_bytes = f.read(space_len)
    times[i] = struct.unpack('<i', leading_bytes[4:8])[0]  # time of measurement
    if wls is None:
        wl_start, wl_end, wl_step = struct.unpack('<HHH', leading_bytes[8:14])
        wls = np.arange(wl_start, wl_end + wl_step, wl_step) * wl_scale
        data_mat = np.empty((nrec, wls.shape[0]), dtype=np.float64)  # create a data matrix for our data

    data_mat[i, :] = np.frombuffer(f.read(8 * data_mat.shape[1]), dtype='<d')  # read the values of fluorescence spectra

data_mat *= scale_fac
times *= scale_time
print(data_mat)
print(data_mat.shape)
print(times)
print(wls)

save_mat2csv(r'C:\Users\Dominik\Documents\Python + JS\File-Converters\test.csv', data_mat, times, wls)
# save_bytes(r'C:\Users\Dominik\Documents\Python + JS\File-Converters\test', b''.join(b_array))

f.close()


# n_records = int((0x4AE10 - 0x4A638) / 8)
# data1 = np.frombuffer(f.read(8 * n_records), dtype='<d')  # read little endian float64 values
# print(data1)
# print(data1.shape)


# FLD1A,Ex=300, Em=430
# < TimeStart > 324 < / TimeStart >
# < TimeEnd > 1800000 < / TimeEnd >
# < Minimum > 0 < / Minimum >
# < Maximum > 0.066898 < / Maximum >
# < Slope > 1E-06 < / Slope >
# < NumberOfValues > 4166 < / NumberOfValues >


# f.seek(0x2D77A2)
# n_records = int((0x2DF9D2 - 0x2D77A2) / 8)
# fld_trace = np.frombuffer(f.read(8 * n_records), dtype='<d')  # read little endian float64 values
# print(fld_trace)
# print(fld_trace.shape)


if __name__ == '__main__':
    a = 0
    # print(str0)
    # print(str1)
    # print(data)
    # print(diffs)
    # plt.plot(data1[::-1])
    # plt.show()
    # print(data1.shape)










