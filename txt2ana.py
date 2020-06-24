#!/usr/bin/python3

# txt2ana.py -- Convert a .txt file from text format to OPTIMUS .ana input file

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
        _data = np.genfromtxt(filename, delimiter='\t', dtype=np.float64)
        times = _data[0, 1:]
        wavelengths = _data[1:, 0]
        D = _data[1:, 1:].T

        fpath_wout_ext = os.path.splitext(filename)[0]
        dir, basename = os.path.split(fpath_wout_ext)
        basename += '.ana'

        delimiter = ' '
        
        timelist = delimiter.join(f'{t}' for t in times)
        wlist = delimiter.join(f'{w}' for w in wavelengths)
        
        buffer = f'%FILENAME={basename}\n%DATATYPE=TAVIS\n%NUMBERSCANS=1\n%TIMESCALE=ps\n%TIMELIST={timelist}\n%WAVELENGTHLIST={wlist}\n%INTENSITYMATRIX=\n'
        
        buffer += '\n'.join(delimiter.join(f"{num}" for num in row) for row in D)

        with open(os.path.join(dir, basename), 'w', encoding='utf8') as f:
            f.write(buffer)
exit(0)
