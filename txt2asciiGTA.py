#!/usr/bin/python3

# txt2asciiGTA.py -- Convert a .txt file from SPECFIT format to Glotaran .ascii file

import sys
import os
import numpy as np

if __name__ == '__main__':
    # Show file open dialog if no input files specified on command line
    if len(sys.argv) > 1:
        filenames = sys.argv[1:]
    else:
        from tkinter import Tk, filedialog
        tk_root = Tk()
        tk_root.withdraw()
        filenames = list(filedialog.askopenfilenames(parent=None, title='Select csv input files...'))
        tk_root.destroy()
        if len(filenames) < 1:
            print('No input files selected. Exiting.')
            exit(1)

    for filename in filenames:
        _data = np.loadtxt(filename, delimiter='\t', skiprows=3, dtype=np.float64)
        times = _data[0, 1:] * 1e-3
        wavelengths = _data[1:, 0]
        D = _data[1:, 1:].T

        new_fname = os.path.splitext(filename)[0]

        delimiter = '\t'

        mat = np.vstack((wavelengths, D))
        buffer = f'Header\nOriginal filename: fname\nTime explicit\nintervalnr {times.shape[0]}\n'
        buffer += delimiter + delimiter.join(f"{num}" for num in times) + '\n'
        buffer += '\n'.join(delimiter.join(f"{num}" for num in row) for row in mat.T)

        with open(new_fname + '-GLOTARAN.ascii', 'w', encoding='utf8') as f:
            f.write(buffer)
exit(0)
