#!/usr/bin/python3

# Converts the Agilent *.UV files to csv, uses 3rd party Aston library https://github.com/bovee/Aston
# works on files form old HPLC
# does not work for new HPLC, files exported as old format (*.D Chemstation) are not read properly for some reason



import sys
import os
import numpy as np

# from https://github.com/bovee/Aston
from aston.tracefile import TraceFile


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

        # from https://github.com/bovee/Aston
        tf = TraceFile(filename)

        chromatogram = tf.data  # https://github.com/bovee/Aston/blob/315871346df72b3e8fcfa9943e8a3519e60299ff/aston/tracefile/__init__.py#L68

        elutionTimes = np.asarray(chromatogram.index)
        wavelengths = np.asarray(chromatogram.columns)
        data = chromatogram.values  # shape: (el. times x wavelengths)

        fpath_wout_ext = os.path.splitext(filename)[0]
        dir, basename = os.path.split(fpath_wout_ext)
        basename += '.csv'

        delimiter = ','

        mat = np.hstack((elutionTimes[:, None], data))
        buffer = f'Elution Time | Wavelength->'
        buffer += delimiter + delimiter.join(f"{num}" for num in wavelengths) + '\n'
        buffer += '\n'.join(delimiter.join(f"{num}" for num in row) for row in mat)
        
        with open(os.path.join(dir, basename), 'w', encoding='utf8') as f:
            f.write(buffer)
        
exit(0)

