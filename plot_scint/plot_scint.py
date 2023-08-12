import numpy as np
from collections import deque
from itertools import islice
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import click
import glob
import os
import re
from multiprocessing import Pool
import fil_spec as fsp

#################
##  SMOOTHING  ##
#################

def moving_average(dat, n=100):
    it = iter(dat)
    d = deque(islice(it, n-1))
    d.appendleft(0)
    s = np.sum(d)
    for elem in it:
        s += elem - d.popleft()
        d.append(elem)
        yield s / float(n)

def moving_average_full(dat, n=100):
    ma = np.array([ ii for ii in moving_average(dat, n) ])
    front_pad = (len(dat) - len(ma)) // 2
    back_pad  = len(dat) - front_pad - len(ma)
    ma = np.hstack( (np.hstack((np.zeros(front_pad), ma)),
                     np.zeros(back_pad)) )
    return ma

def rebin(x, nbins=1024, strip_pad=False):
    avgkern = int(len(x) / nbins)
    avgx = moving_average_full(x, n=avgkern)
    tt0 = np.linspace(0, 1, len(avgx))
    tt1 = np.linspace(0, 1, nbins)
    nx = np.interp(tt1, tt0, avgx)
    if strip_pad:
        if nx[0] == 0:
            nx = nx[1:]
        if nx[-1] == 0:
            nx = nx[:-1]
    return nx

def rebin_spec(freqs, spec, dsfac=1, strip_pad=True):
    N = len(freqs)
    nbins = int( N / dsfac )
    afreqs = rebin(freqs, nbins=nbins, strip_pad=strip_pad)
    aspec  = rebin(spec, nbins=nbins, strip_pad=strip_pad)
    return afreqs, aspec

def rebin2(x, avgkern=1):
    nbins = int(len(x) / avgkern)
    avgx = moving_average_full(x, n=avgkern)
    tt0 = np.linspace(0, 1, len(avgx))
    tt1 = np.linspace(0, 1, nbins)
    nx = np.interp(tt1, tt0, avgx)
    return nx


##########
## ACFS ##
##########

def calc_acf(x):
    a = np.correlate(x, x, mode='full')
    return a[len(a)//2:]

def calc_ccf(x, z):
    a = np.correlate(x, z, mode='full')
    N = len(z)
    lags = np.arange(-1 * N + 1, N, 1)
    return lags, a

def fft_conv(a, b):
    cc = np.fft.irfft( np.fft.rfft(a) * np.fft.rfft(b) )
    cc /= np.sum(b)
    return cc

def calc_norm_acf(x, y, normbin=-1):
    
    a = calc_acf(y)
    lags = np.abs(x - x[0])
    if normbin > -1:
        a = a / a[int(normbin)]
    return lags, a

##############
## SPEC ACF ##
##############

def spec_acf(freqs, spec, normbin=1):
    """
    calculate spectrum acf 
    """
    m = np.mean(spec)
    lags, acf = calc_norm_acf(freqs, spec-m, normbin=normbin)
    return lags, acf


def pulse_acf(tt, prof, normbin=1):
    """
    calculate pulse acf 
    """
    m = np.mean(prof)
    lags, acf = calc_norm_acf(tt, prof-m, normbin=normbin)
    return lags, acf


def simple_scint_bw(freqs, spec, favg=0.3):
    """
    Calc scint bw as 1/e point of ACF

    subtract a smoothed version of the spectrum 
    that is averaged over favg * bw 
    """
    # Calc smoothed spec
    navg = int( favg * len(spec) )  
    if navg < 1:
        navg = 1
    else: pass
    spec_avg = moving_average_full(spec, navg)
    
    # Calc ACF and normalize so ACF(0) = 1
    lags, acf = spec_acf(freqs, spec-spec_avg, normbin=0)

    # Find first instance where ACF(x) <= 1/e
    xx = np.where( acf <= 1/np.e )[0][0]

    sbw = lags[xx] 

    return lags, acf, sbw, spec_avg


def simple_pulse_acf(tt, prof):
    """
    Calc char width as 1/e point of ACF
    """
    # Calc ACF and normalize so ACF(0) = 1
    lags, acf = pulse_acf(tt, prof, normbin=0)

    # Find first instance where ACF(x) <= 1/e
    xx = np.where( acf <= 1/np.e )[0][0]

    w = lags[xx] 

    return lags, acf, w


def scint_plot(freqs, spec, out_file):
    """
    Make a plot of the spectrum + ACF
    """
    fig = plt.figure(figsize=(6,6))
    
    ax_s = fig.add_subplot(311)
    ax_a = fig.add_subplot(312)
    ax_z = fig.add_subplot(313)

    lags, acf, sbw, savg = simple_scint_bw(freqs, spec)

    fs = 12

    # Spec plot
    ax_s.plot(freqs, spec)  
    smax = np.max(spec)
    fmid = np.mean(freqs)
    ax_s.errorbar(fmid, 0.9 * smax, xerr=sbw/2., 
                  marker='', capsize=3, capthick=1.5, 
                  color='r')
    ax_s.plot(freqs, savg, c='k', ls='--')

    ax_s.set_xlim(np.min(freqs), np.max(freqs))
    ax_s.set_xlabel("Frequency (MHz)", fontsize=fs)
    ax_s.set_ylabel("Flux Density (arb)", fontsize=fs)
    

    # ACF -- Full
    ax_a.plot(lags, acf, marker='o')
    ax_a.set_xlim(0, lags[-1])
    #ax_a.axvline(x=sbw, ls='--', c='k')
    ax_a.set_xlabel("Frequency Lag (MHz)", fontsize=fs)
    ax_a.set_ylabel("ACF", fontsize=fs)

    #sbw_str = r"$\Delta \nu_{\rm s} = %.1f\, {\rm MHz}$" %sbw
    #ax_a.text(0.9, 0.8, sbw_str, ha='right', va='top', 
    #          transform=ax_a.transAxes, fontsize=16)

    # ACF -- zoom
    dl = np.diff(lags)[0]
    xmax = max(10 * dl, min(sbw * 5, lags[-1]))
    ax_z.plot(lags, acf, marker='o')
    ax_z.set_xlim(0, xmax)
    ax_z.axvline(x=sbw, ls='--', c='k')
    ax_z.set_xlabel("Frequency Lag (MHz)", fontsize=fs)
    ax_z.set_ylabel("ACF", fontsize=fs)
   
    vstr, ustr = get_sbw_string(sbw) 
    lab_str = r"$\Delta \nu_{\rm s} = %s\, {\rm %s}$" %(vstr, ustr)
    ax_z.text(0.9, 0.8, lab_str, ha='right', va='top', 
              transform=ax_z.transAxes, fontsize=16)

    plt.savefig(out_file)

    return sbw


def get_sbw_string(sbw):
    """
    format sbw string 
    """
    if sbw <= 0:
        vstr = "%.1f" %sbw
        ustr = "MHz"
        return vstr, ustr
    else: pass
    
    log_s = np.log10(sbw)
    if log_s >= 3:
        vstr =  "%.1f" %(sbw/1e3)
        ustr = "GHz"

    elif 0 <= log_s < 3:
        vstr = "%.1f" %sbw
        ustr = "MHz"

    elif -3 <= log_s < 0:
        vstr = "%.1f" %(sbw * 1e3)
        ustr = "kHz"

    else:
        vstr = "%.1f" %(sbw * 1e6)
        ustr = "Hz"

    return vstr, ustr


###############################
## Processing multiple files ##
###############################

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


def generate_plots(offsets, sp_dt, bin_widths, infile,  
                   png_dir, dm, tavg, tdur):
    """
    Creates *.png plots of spectra from fil files
    Takes in offset times to find where to plot
    """
    # Get corresponding offset value by checking which pulse num
    filename = os.path.basename(infile)

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

    pulse_width = bin_width * sp_dt * 1e-6

    basename = os.path.splitext(filename)[0]
    # plot using sp_spec

    freqs, spec = fsp.get_spec(infile, dm, tavg=int(tavg), tp=offset, 
                               tdur=tdur, pulse_width=pulse_width, 
                               bpass=False)

    #np.save('freqs.npy', freqs)
    #np.save('spec.npy', spec)
    
    # extract basename
    basename = os.path.basename(infile)
    basename = os.path.splitext(basename)[0]
    outfile = '%s/%s_scint.png' % (png_dir, basename)

    sbw = scint_plot(freqs, spec, outfile)

    return
     

@click.command()
@click.option("--off_file", type=str,
              help="Full path of offset file", required=True)
@click.option("--fil_dir", type=str,
              help="Dir of *fil files", required=True)
@click.option("--comb_base", type=str,
              help="Combined filterbank basename", required=True)
@click.option("--png_dir", type=str,
              help="Dir of *png file", required=True)
@click.option("--dm", type=float,
              help="Disperion Measure", required=True)
@click.option("--tavg", type=float,
              help="Average time", required=True)
@click.option("--tdur", type=float,
              help="Duration of pulse", required=True)
@click.option("--sp_dt", type=float,
              help="Time res of SP file (us)", required=True)
@click.option("--nproc", type=int, default=1,
              help="Number of processes for multiprocessing", 
              required=False)
def plot_spec(off_file, sp_dt, fil_dir, comb_base, 
              png_dir, dm, tavg, tdur, nproc=1):
    """
    Creates plot of spectrum and ACF for scintillation study
    """
    offsets, snrs, bin_widths, tcands = read_data(off_file)
    fil_files = glob.glob('%s/%s*fil' % (fil_dir, comb_base))

    fil_files.sort()

    if len(fil_files) == 0:
        print("No fil files found!")
        print("fil_dir: %s" %fil_dir)
        print("fil_base: %s" %comb_base)
        return 0
    else: pass

    #for ii, ff in enumerate(fil_files):
    #    generate_plots(offsets, bin_widths, ff, png_dir, 
    #                   dm, tavg, tdur)

    with Pool(nproc) as pool:
        args_list = [(offsets, sp_dt, bin_widths, fil_file,
                      png_dir, dm, tavg, tdur) for fil_file in fil_files]
        pool.starmap(generate_plots, args_list)


debug = 0

if __name__ == "__main__":
    if not debug:
        plot_spec()
