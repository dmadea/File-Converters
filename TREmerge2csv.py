#!/usr/bin/python3

# Converts the transient emission csv file to matrix-csv file
# file in the format of e.g. 100acc-0.45exp-air.csv


# import sys
import os
import numpy as np
import re
import argparse
from glob import glob

from hplc_converter_base import save_mat2csv, fi


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

    for i in range(ns):
        exp_time_from_data = data[i * nw, -2]   # exposure time
        _exp_time = exp_time if np.isnan(exp_time_from_data) else exp_time_from_data   # if exposure time is not present in the datafile, extracted value from filename is used
        times[i] = data[i * nw, -1]
        D[i, :] = data[i*nw:(i+1)*nw, 0] / (_exp_time * acc)   # divide by exposure time and by number of accumulations used
        # D[i, :] = data[i*nw:(i+1)*nw, 0] / acc   # divide by exposure time and by number of accumulations used

    print(f"File: {fname},\nExposure time: {exp_time if np.isnan(exp_time_from_data) else "varying"}, Number of accumulations: {acc}\n")

    return Dataset(times, wavelengths, D)

def stack_datasets(datasets: list, t0, w0, w1, adjust_area, k, bd_corr, bd_w0, bd_w1):  # , adjust_next=False


    mul_factor = 1
    # x = datasets[0].w
    idx = 1   #  in default it will remove the first scan from all datasets (except the first one)
    # t0 = datasets[0].t[0] if t0 is None else t0 # time zero for the first dataset
    w = datasets[0].w

    # subtract the time zero from all of the datasets time points
    if t0 is not None:
        first_t = datasets[0].t[0]
        for d in datasets:
           d.t -= first_t + t0

    # baseline drift correct datasets
    if bd_corr:
        idx1 = fi(w, bd_w0) if bd_w0 is not None else 0
        idx2 = fi(w, bd_w1) + 1 if bd_w1 is not None else w.shape[0]

        for d in datasets:
            D_selection = d.D[:, idx1:idx2]
            d.D -= D_selection.mean(axis=1, keepdims=True)
    
    for i in range(1, len(datasets)):

        if adjust_area:
        #  make an average of last k spectra
            lastsp = np.trapz(datasets[i - 1].D[-k:, :].mean(axis=0), w)

            sp = np.trapz(datasets[i].D[0, :], w)  # only first spectrum
            mul_factor = lastsp / sp

        datasets[i].D = mul_factor * datasets[i].D[idx:, :]
        datasets[i].t = datasets[i].t[idx:]
        print(f"Dataset {i} mul factor {mul_factor}")

    new_t = np.hstack(list(map(lambda d: d.t, datasets)))
    new_D = np.vstack(list(map(lambda d: d.D, datasets)))

    idx1 = fi(w, w0) if w0 is not None else 0
    idx2 = fi(w, w1) + 1 if w1 is not None else w.shape[0]

    Dcrop = new_D[:, idx1:idx2]
    wcrop = w[idx1:idx2]
        
    return Dataset(new_t, wcrop, Dcrop)


def process_filepaths(fpaths, t0, w0, w1, adjust_area, k, proc1by1, bd_corr, bd_w0, bd_w1):
    print(f"\nProcessing {len(fpaths)} files:")
    for f in fpaths:
        print(f"- {os.path.basename(f)}")
    print(f"\nSettings:")
    print(f"- Time zero adjustment: {t0 if t0 is not None else 'None'}")
    print(f"- Wavelength range: {w0 if w0 is not None else 'start'} to {w1 if w1 is not None else 'end'}")
    print(f"- Area adjustment: {adjust_area}")
    print(f"- Baseline drift correction: {bd_corr}")
    if bd_corr:
        print(f"  - Baseline drift wavelength range: {bd_w0} to {bd_w1}")
    print(f"- Process one by one: {proc1by1}\n")

    # sort the filepaths according to its filename
    fpaths = sorted(fpaths, key=lambda entry: os.path.split(entry)[1])

    #  load datasets

    opt = "_area_ajusted" if adjust_area else ""

    if proc1by1:
        for path in fpaths:
            print(f"\nProcessing file: {os.path.basename(path)}")
            dataset = load_TRE(path)
            mD = stack_datasets([dataset], t0, w0, w1, adjust_area, k, bd_corr, bd_w0, bd_w1)
            _dir, fname = os.path.split(path)   # get dir and filename
            fname, _ = os.path.splitext(fname)
            output_path = os.path.join(_dir, f'{fname}_proc{opt}.csv')
            save_mat2csv(output_path, mD.D, mD.t, mD.w, unit='ns')
            print(f"Saved processed file to: {output_path}")

    else:
        print("\nMerging all datasets...")
        datasets = list(map(lambda fpath: load_TRE(fpath), fpaths))

        mD = stack_datasets(datasets, t0, w0, w1, adjust_area, k, bd_corr, bd_w0, bd_w1)

        _dir, fname = os.path.split(fpaths[0])   # get dir and filename
        fname, _ = os.path.splitext(fname)
        output_path = os.path.join(_dir, f'{fname}-{len(datasets)}_merge{opt}.csv')
        save_mat2csv(output_path, mD.D, mD.t, mD.w, unit='ns')
        print(f"Saved merged file to: {output_path}")


if __name__ == '__main__':
    # Show file open dialog if no input files specified on command line

    parser = argparse.ArgumentParser()

    parser.add_argument("--t0", nargs="?", default=None, type=float,
                        help="Set time zero for the first dataset if proc_1by1==False or for all datasets if proc_1by1==True. If None (default), no change of time zero will be made.")

    parser.add_argument("--w0", nargs="?", default=None, type=float,
                        help="Start wavelength to crop the data.")
    
    parser.add_argument("--w1", nargs="?", default=None, type=float,
                        help="End wavelength to crop the data.")

    parser.add_argument("--adjust_area", action="store_true",
                        help="If True, connects the datasets so that the area of the average of last k-th spectra of n-th dataset is the same as area of first spectrum of n+1-th dataset.")
    
    parser.add_argument("--k", nargs="?", default=5, type=int,
                        help="Number of spectra to average (only work if adjust_area==True).")
    
    parser.add_argument("--proc1by1", action="store_true",
                        help="If True, processes and saved the datasets individually with added '_proc' sufix.")
    
    parser.add_argument("--bd_corr", action="store_true",
                        help="If True, baseline drift correction based on noise region will be performed for each dataset.")
    
    parser.add_argument("--bd_w0", default=100, type=float,
                        help="Start wavelength of the noise region for baseline drift correction.")
    
    parser.add_argument("--bd_w1", default=300, type=float,
                        help="End wavelength of the noise region for baseline drift correction.")

    parser.add_argument('files', nargs=argparse.ZERO_OR_MORE)

    args, _ = parser.parse_known_args()

    filepaths = []
    for fname in args.files:
        filepaths += glob(fname)

    if len(filepaths) == 0:
        from tkinter import Tk, filedialog

        tk_root = Tk()
        tk_root.withdraw()
        filepaths = list(filedialog.askopenfilenames(parent=None, title='Select csv input files from TRE data.'))
        tk_root.destroy()
        if len(filepaths) < 1:
            print('No input files selected. Exiting.')
            exit(1)


    process_filepaths(filepaths, args.t0, args.w0, args.w1, args.adjust_area, args.k, args.proc1by1, args.bd_corr, args.bd_w0, args.bd_w1)