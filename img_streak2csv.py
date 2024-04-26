#!/usr/bin/python3

# Converts the transient emission csv file to matrix-csv file
import sys
import os
import pandas
from streakimage import StreakImage

def process_filepaths(fpaths):
    for path in fpaths:
        img = StreakImage(path)
        _dir, fname = os.path.split(path)   # get dir and filename
        fname, _ = os.path.splitext(fname)
        img.data.iloc[:, ::-1].to_csv(os.path.join(_dir, f'{fname}.csv'))  # reverse the columns and export to csv

if __name__ == '__main__':
    # Show file open dialog if no input files specified on command line
    if len(sys.argv) > 1:
        filepaths = sys.argv[1:]
    else:
        from tkinter import Tk, filedialog

        tk_root = Tk()
        tk_root.withdraw()
        filepaths = list(filedialog.askopenfilenames(parent=None, title='Select csv input files from TRE data.'))
        tk_root.destroy()
        if len(filepaths) < 1:
            print('No input files selected. Exiting.')
            exit(1)

    process_filepaths(filepaths)