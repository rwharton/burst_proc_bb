import numpy as np
import os 
import glob
import time
import subprocess
import click
import re
from datetime import datetime, timedelta
from astropy.time import Time, TimeDelta
import sp_spec as sp
import multiprocessing as mp
import matplotlib.pyplot as plt


def read_data(offset_file):
    """
    Read in offset, SNR
    """
    dd = np.loadtxt(offset_file)
    offsets = dd[:, 1]
    widths  = dd[:, 2].astype('int')
    snrs    = dd[:, 3]
    times   = dd[:, 4]  
    return offsets, snrs, widths, times


def generate_plots(offsets, sp_dt, snrs, bin_widths, tcands, 
                   infile, comb_base, png_dir, dm, plt_dt, 
                   plt_chans, tdur):
    """
    Creates *.png plots of filterbank files
    Takes in offset times to find where to plot
    """
    # Get corresponding offset value by checking which pulse num
    filename = os.path.basename(infile)

    # Use regex to find the pattern. We are looking for p 
    # followed by 4 digits
    match = re.search(r'p(\d{4})', filename)
    index = -1
    if match:
        # Convert the digit part to an integer
        index = int(match.group(1))
    else:
        print("Pattern not found in filename.")

    offset = offsets[index]
    snr = snrs[index]
    bin_width = bin_widths[index]
    tcand = tcands[index]

    dt = sp_dt * 1e-6

    pulse_width = bin_width * dt

    basename = os.path.splitext(filename)[0]
    # plot using sp_spec

    # no need to bandpass
    bpass = True
    if comb_base in basename:
        bpass = False
        
    freqs, spec = sp.make_plot('%s' % infile, dm, index, snr, 
                               plt_dt=plt_dt, plt_chans=plt_chans, 
                               tp=offset, tdur=tdur, 
                               outbins=1024, tcand=tcand, 
                               outfile="%s/%s.png" % (png_dir, basename), 
                               pulse_width=pulse_width, bpass=bpass)
    
    return 


@click.command()
@click.option("--off_file", type=str, 
              help="Full path of offset file", required=True)
@click.option("--sp_dt", type=float, 
              help="Time res of SP file (us)", required=True)
@click.option("--fil_dir", type=str, 
              help="Dir of *fil files", required=True)
@click.option("--fil_base", type=str,
              help="basename of *fil files", required=True)
@click.option("--comb_base", type=str,
              help="Combined filterbank basename", required=True)
@click.option("--png_dir", type=str, 
              help="Dir of *png file", required=True)
@click.option("--dm", type=float, 
              help="Disperion Measure", required=True)
@click.option("--plt_dt", type=float, 
              help="Plot time res (ms)", required=True)
@click.option("--plt_chans", type=int, 
              help="Avg to this number of channels", 
              required=True)
@click.option("--tdur", type=float, 
              help="Data length for plot", required=True)
@click.option("--nproc", type=int, default=1, 
              help="Number of processes for multiprocessing", 
              required=False)
def plotfil(off_file, sp_dt, fil_dir, fil_base, comb_base, 
            png_dir, dm, plt_dt, plt_chans, tdur, nproc=1):
    """
    Creates plot of filterbank files.

    Takes in filterbank directory and plots burst given an 
    offset value for each file.
    """
    offsets, snrs, bin_widths, tcands = read_data(off_file)
    fil_files = glob.glob('%s/%s*fil' % (fil_dir, fil_base))

    if len(fil_files) == 0:
        print("No fil files found!")
        print("fil_dir: %s" %fil_dir)
        print("fil_base: %s" %fil_base)
        return 0
    else: pass

    #for infile in fil_files:
    #    generate_plots(offsets, sp_dt, snrs, bin_widths, infile, 
    #                   comb_base, png_dir, dm, tavg, tdur)

    with mp.Pool(nproc) as pool:
        args_list = [(offsets, sp_dt, snrs, bin_widths, tcands, 
                      infile, comb_base, png_dir, dm, plt_dt, 
                      plt_chans, tdur) for infile in fil_files]
        pool.starmap(generate_plots, args_list)
    return

debug = 0

if __name__ == "__main__":
    if not debug: 
        plotfil()


