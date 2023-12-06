#!/usr/bin/python3

# Converts the transient emission csv file to matrix-csv file
import sys
import os
import numpy as np

from hplc_converter_base import save_mat2csv

def load_TRE(filename, nw=1024):

    data = np.genfromtxt(filename, delimiter=',', dtype=np.float64, usecols=(4, 5, 10))
    ns = int(data.shape[0] / nw)
    
    times = np.empty(ns, dtype=float)
    wavelengths = data[:nw, 1]
    D = np.empty((ns, nw), dtype=float)
    
    for i in range(ns):
        times[i] = data[i * nw, -1]
        D[i, :] = data[i*nw:(i+1)*nw, 0]

    return times, wavelengths, D

def process_filepath(fpath):

    _dir, fname = os.path.split(fpath)   # get dir and filename
    fname, _ = os.path.splitext(fname)  # get filename without extension

    t, w, D = load_TRE(fpath)
    save_mat2csv(os.path.join(_dir, f'{fname}-export.csv'), D, t, w, unit='ns')



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

