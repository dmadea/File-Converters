#!/usr/bin/python3

# Converts the transient emission csv file to matrix-csv file
import sys
import os
import numpy as np
import re

from hplc_converter_base import save_mat2csv

racc = re.compile(r"(\d+)acc")
rexposure = re.compile(r"(\d+\.?\d*)exp")


class Dataset:
    def __init__(self, t, w, D):
        self.t = t
        self.w = w
        self.D = D

def load_TRE(filename, nw=1024):
    
    _dir, fname = os.path.split(filename)   # get dir and filename
    fname, _ = os.path.splitext(fname)  # get filename without extension

    data = np.genfromtxt(filename, delimiter=',', dtype=np.float64, usecols=(4, 5, 9, 10))
    ns = int(data.shape[0] / nw)
    
    times = np.empty(ns, dtype=float)
    wavelengths = data[:nw, 1]
    D = np.empty((ns, nw), dtype=float)

    macc = racc.search(fname)
    mexp = rexposure.search(fname)

    acc = 1  # number of accumulations
    exp_time = 1
    if macc is not None:
        acc = int(macc.group(1))

    if mexp is not None:
        exp_time = float(mexp.group(1))
    
    print(f"File: {fname},\nExposure time: {exp_time}, Number of accumulations: {acc}\n")

    for i in range(ns):
        exp_time_from_data = data[i * nw, -2]   # exposure time
        exp_time = exp_time if np.isnan(exp_time_from_data) else exp_time_from_data   # if exposure time is not present in the datafile, 0.45 value is used
        times[i] = data[i * nw, -1]
        D[i, :] = data[i*nw:(i+1)*nw, 0] / (exp_time * acc)   # divide by exposure time and by number of accumulations used
        # D[i, :] = data[i*nw:(i+1)*nw, 0] / acc   # divide by exposure time and by number of accumulations used


    return Dataset(times, wavelengths, D)

def stack_datasets(*datasets):  # , adjust_next=False

    mul_factor = 1
    # x = datasets[0].w
    idx = 1   #  in default it will remove the first scan from all datasets (except the first one)
    #t0 = datasets[0].t[0]  # time zero for the first dataset

    # subtract the time zero from all of the datasets time points
    #for d in datasets:
    #    d.t -= t0
    
    for i in range(1, len(datasets)):
        datasets[i].t = datasets[i].t[idx:]

        # lastsp = np.trapz(datasets[i - 1].D[-1, :], x)
        # sp = np.trapz(datasets[i].D[0, :])
        # mul_factor *= lastsp / sp
        datasets[i].D = mul_factor * datasets[i].D[idx:, :]

    new_t = np.hstack(list(map(lambda d: d.t, datasets)))
    new_D = np.vstack(list(map(lambda d: d.D, datasets)))
        
    return Dataset(new_t, datasets[0].w, new_D)

def process_filepaths(fpaths):

    # sort the filepaths according to its filename
    fpaths = sorted(fpaths, key=lambda entry: os.path.split(entry)[1])

    for path in fpaths:
        dataset = load_TRE(path)
        mD = stack_datasets(dataset)
        _dir, fname = os.path.split(path)   # get dir and filename
        fname, _ = os.path.splitext(fname)
        save_mat2csv(os.path.join(_dir, f'{fname}_proc.csv'), mD.D, mD.t, mD.w, unit='ns')


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