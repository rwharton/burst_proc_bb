import sp_spec as sp
import numpy as np
import click
import multiprocessing as mp
import pandas as pd
import os
import re
import matplotlib.pyplot as plt
import glob


def read_offset(fil_file, offset_file):
    """
    Read in offset
    """
    df = pd.read_csv(offset_file)
    offsets = df["Offset Times"]
    bin_widths = df["Bin Width"]
        # Get corresponding offset value by checking which pulse num
    filename = os.path.basename(fil_file)

    # Use regex to find the pattern. We are looking for p followed by 4 digits
    match = re.search(r'p(\d{4})', filename)
    index = -1
    if match:
        # Convert the digit part to an integer
        index = int(match.group(1))
    else:
        print("Pattern not found in filename.")
    offset = offsets[index]
    bin_width = bin_widths[index]
    return offset, bin_width


def calc_snr(fil_file, dm, pulse_bins, tp, tdur, bin_width = 1.024*10**(-5), sample_time = 10**(-6), n_bins=5):
    """
    Extra time series data from filterbank and calculate snr
    pulse_bins: Width of pulse in bins
    bin_width: Length of bin time in seconds
    n_bins: new width of pulse in bins
    """
    # number of bins to avg over
    tavg = ((pulse_bins * bin_width) / n_bins) / sample_time
    bpass = True
    if "bnd" in fil_file:
        bpass = False
    tt, freqs, dat = sp.get_snippet_data(fil_file, dm, 
                             favg=0, tavg=int(tavg), bpass=bpass,
                             tp=tp, tdur=tdur)
    
    # get ts from dat
    ts = np.mean(dat, axis=1)
    Nthird = len(ts) // 3
    avg = np.mean(ts[:Nthird])
    sig = np.std(ts[:Nthird])
    ts = (ts-avg)/sig

    # get out of signal range
    shift = 100

    signal = np.max(ts)
    noise = np.std(ts[:np.argmax(ts)-shift])

    return signal/noise, dm

def write_dms(out_dir, dm_dir, basename, dms):
    """
    Write results to out file and generates png
    """
    df = pd.DataFrame(dms, columns = ['SNR', 'DM'])
    # reorder columns to put DM first

    df = df[['DM', 'SNR']]
    df.to_csv('%s%s%s.dm' % (out_dir, dm_dir, basename), index=False, sep='\t')

    # create graph
    df = df.sort_values('DM')
    plt.figure()
    plt.plot(df.iloc[:, 0], df.iloc[:, 1]) 
    plt.xlabel('DM') 
    plt.ylabel('SNR') 
    plt.savefig("%s%s%s.png" % (out_dir, dm_dir, basename))




@click.command()
@click.option("--out_dir", type=str,
              help="Path to out dir", required=True)
@click.option("--dm_dir", type=str,
              help="Path to dm dir", required=True)
@click.option("--fil_dir", type=str, 
              help="Path to filterbank dir", required=True)
@click.option("--offsets", type=str,
              help="Path to offsets file", required=True)
@click.option("--dm_lo", type=float, 
              help="Lower bound of dm to search for", required=True)
@click.option("--dm_hi", type=float,
              help="Upper bound of dm to search for", required=True)
@click.option("--dm_step", type=float, 
              help="Step value of dm", required=True)

def dm_opt(out_dir, dm_dir, fil_dir, offsets, dm_lo, dm_hi, dm_step):
    """
    Optimizes dm by finding maximum SNR
    """
    fil_files = glob.glob("%s%s*.fil" % (out_dir, fil_dir))
    print(fil_files)
    for fil_file in fil_files:
        offset, bin_width = read_offset(fil_file, out_dir+offsets)
        print(offset, bin_width)
        inputs = [(fil_file, round(dm, 2), bin_width, round(offset, 2), 0.2) for dm in np.arange(dm_lo, dm_hi, dm_step)]


        with mp.Pool() as pool:
            results = pool.starmap(calc_snr, inputs)

        max_snr, max_dm = max(results)
        print("Maximum SNR: %f at DM: %f" %(max_snr, max_dm))
        
        base_name = os.path.basename(fil_file)
        # Get the file name without the extension.
        base_name= os.path.splitext(base_name)[0]

        # write in sorted out
        write_dms(out_dir, dm_dir, base_name, sorted(results, reverse=True))



debug = 0
if __name__ == "__main__":
    if not debug:
        dm_opt()
