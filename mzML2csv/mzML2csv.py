#!/usr/bin/python3

# Converts the mzML files to csv, uses 3rd party pyteomics library https://github.com/levitsky/pyteomics 
# for parsing the mzML file, the only parameter that can be changed here in this file is resolution of MS data
# default value is 1, change is if necessary

# For conversion of mzML files from *.wiff and *wiff.scan files, use online converter https://gnps-quickstart.ucsd.edu/conversion

import sys
import os
import numpy as np
from pyteomics.mzml import read, iterfind

MZ_RESOLUTIONs = [0.2, 0.5, 1.0]  # change the value here


def save_mat2csv(fname, matrix, times=None, wls=None, delimiter=','):
    times = np.arange(0, matrix.shape[0]) if times is None else times
    wls = np.arange(0, matrix.shape[1]) if wls is None else wls

    mat = np.hstack((times[:, None], matrix))
    buffer = f'MS data - Time | m/z->'
    buffer += delimiter + delimiter.join(f"{num}" for num in wls) + '\n'
    buffer += '\n'.join(delimiter.join(f"{num}" for num in row) for row in mat)

    with open(fname, 'w', encoding='utf8') as f:
        f.write(buffer)

def parse_mzml_file(filepath, mz_resolution=1):

    mzml_file = read(filepath)
    spectrumList_dict = iterfind(filepath, 'indexedmzML/mzML/run/spectrumList', read_schema=True, recursive=False).__next__()
    n_spectra = int(spectrumList_dict['count'])  # number of all spectra

    mzmin, mzmax = None, None
    mz_array = None
    mat = None
    times = None

    for i, sp in enumerate(mzml_file):

        if mat is None:
            mzmin = float(sp['scanList']['scan'][0]['scanWindowList']['scanWindow'][0]['scan window lower limit'])
            mzmax = float(sp['scanList']['scan'][0]['scanWindowList']['scanWindow'][0]['scan window upper limit'])

            mzmax = mzmax + mz_resolution - (mzmax - mzmin) % mz_resolution  # recalculate the maximum m/z value

            mz_array = np.linspace(mzmin, mzmax, int((mzmax - mzmin) / mz_resolution) + 1) # make sure to have evenly spaced integerers
            mat = np.zeros((n_spectra, mz_array.shape[0]))
            times = np.zeros(n_spectra)
            diff = mzmax - mzmin

        indexes = ((sp['m/z array'] - mzmin) * (mz_array.shape[0] - 1) / diff).astype(int)

        # find duplicated and integrate (just sum) them 
        indexes, indices, counts = np.unique(indexes, return_index=True, return_counts=True)
        intensities = sp['intensity array']
        integrated_intensities = intensities[indices]

        for j in range(indexes.shape[0]):
            if counts[j] < 2:
                continue

            idx = indices[j]
            integrated_intensities[j] = intensities[idx:idx + counts[j]].sum()

        try:
            times[i] = float(sp['scanList']['scan'][0]['scan start time'])
        except KeyError:
            pass
        mat[i, indexes] = integrated_intensities

    return mat, times, mz_array


if __name__ == '__main__':
    # Show file open dialog if no input files specified on command line
    if len(sys.argv) > 1:
        filenames = sys.argv[1:]
    else:
        from tkinter import Tk, filedialog
        tk_root = Tk()
        tk_root.withdraw()
        filenames = list(filedialog.askopenfilenames(parent=None, title='Select mzML input files...'))
        tk_root.destroy()
        if len(filenames) < 1:
            print('No input files selected. Exiting.')
            exit(1)

    for filename in filenames:

        fpath_wout_ext = os.path.splitext(filename)[0]
        dir, basename = os.path.split(fpath_wout_ext)

        for res in MZ_RESOLUTIONs:
            new_filename = f'{basename}-{res}res.csv'
            matrix, times, mz_array = parse_mzml_file(filename, res)
            save_mat2csv(os.path.join(dir, new_filename), matrix, times, mz_array)
        
exit(0)

