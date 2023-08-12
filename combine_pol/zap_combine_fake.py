import numpy as np
import matplotlib.pyplot as plt
import your
from your.formats.filwriter import make_sigproc_object
from scipy import signal
import time
from argparse import ArgumentParser
import os

#######################
###  GET FILE INFO  ###
#######################

def get_fil_info(data_file):
    """
    Get obs info
    """
    yr = your.Your(data_file)
    foff = yr.your_header.foff
    fch1 = yr.your_header.fch1
    dt   = yr.your_header.tsamp
    nchans = yr.your_header.nchans
    tstart = yr.your_header.tstart
    nspec = yr.your_header.nspectra
    return nchans, fch1, foff, dt, tstart, nspec


def check_files(infile1, infile2):
    """
    Check to make sure files are consistent 
    with each other and can be combined
    """
    infoA = get_fil_info(infile1)
    infoB = get_fil_info(infile2)

    nchansA, fch1A, foffA, dtA, tstartA, nspecA = infoA
    nchansB, fch1B, foffB, dtB, tstartB, nspecB = infoB

    err_val = 0

    if nchansA != nchansB:
        print("nchan mismatch: %d, %d" %(nchansA, nchansB))
        err_val += 1
    
    if fch1A != fch1B:
        print("fch1 mismatch: %.2f, %.2f" %(fch1A, fch1B))
        err_val += 1

    if foffA != foffB:
        print("foff mismatch: %.2f, %.2f" %(foffA, foffB))
        err_val += 1
    
    if dtA != dtB:
        print("dt mismatch: %.2f us, %.2f us" %(dtA*1e-6, dtB*1e6))
        err_val += 1

    if nspecA != nspecB:
        print("nspec mismatch: %d, %d" %(nspecA, nspecB))
        err_val += 1
    
    if tstartA != tstartB:
        print("tstart mismatch: %.10f, %.10f" %(tstartA, tstartB))
        err_val += 1

    return err_val

    
###################
###  GET DATA   ###
###################

def get_data_real(infile):
    """
    Read in data
    """
    # get info
    info = get_fil_info(infile)
    nchans, fch1, foff, dt, tstart, nspec = info

    # get times and freqs
    tt = np.arange(nspec) * dt
    freqs = np.arange(nchans) * foff + fch1 

    # get data 
    yr = your.Your(infile)
    dat = yr.get_data(0, nspec)
    #dat = np.random.normal(loc=5, scale=1, size=(37494, 8000))

    return tt, freqs, dat


def get_data(infile):
    """
    "Read" in "data"
    """
    # get info
    tt = np.arange(37494) * 8e-5
    freqs = 2299.99375 + np.arange(8000) * -0.01245
    dat = np.random.normal(loc=5, scale=1, size=(37494, 8000))

    return tt, freqs, dat


####################
###   BP Calcs   ###
####################

def moving_median(data, window):
    """
    Calculate running median and stdev
    """
    startIdxOffset = np.floor(np.divide(window, 2.0))
    endIdxOffset = np.ceil(np.divide(window, 2.0))

    startIndex = startIdxOffset
    endIndex = len(data) - 1 - endIdxOffset

    halfWindow = 0.0

    if (np.mod(window, 2.0) == 0):
        halfWindow = int(np.divide(window, 2.0))
    else:
        halfWindow = int(np.divide(window - 1.0, 2.0))

    mov_median = np.zeros(len(data))
    mov_std = np.zeros(len(data))

    startIndex = int(startIndex)
    endIndex = int(endIndex)

    # Calculate the moving median and std. dev. associated
    # with each interval.
    for i in np.arange(startIndex, endIndex + 1, 1):
        istart = int(i - halfWindow)
        istop  = int(i + halfWindow + 1)
        medianValue = np.median(data[istart : istop])
        stdValue = np.std(data[istart : istop], ddof=1)

        mov_median[i] = medianValue
        mov_std[i] = stdValue

    # Set the values at the end points.
    for i in np.arange(0, startIndex, 1):
        mov_median[i] = mov_median[startIndex]
        mov_std[i] = mov_std[startIndex]

    for i in np.arange(endIndex + 1, len(data), 1):
        mov_median[i] = mov_median[endIndex]
        mov_std[i] = mov_std[endIndex]

    return mov_median, mov_std


def plot_bp(freqs, bp, mask_chans, diff_thresh=None,
            val_thresh=None, outfile=None):
    """
    Plot bandpass with masked chans indicated
    """
    chans = np.arange(0, len(freqs), 1)
    good_chans = np.setdiff1d(chans, mask_chans)

    # if outputting file, turn off interactive mode
    if outfile is not None:
        plt.ioff()
    else:
        plt.ion()

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.plot(freqs[good_chans], bp[good_chans], 'k.')
    ax.plot(freqs[mask_chans], bp[mask_chans], ls='',
            marker='o', mec='r', mfc='none')

    if val_thresh is not None:
        ax.axhline(y=val_thresh, ls='--', c='g', alpha=0.5)
    else: pass

    ax.set_xlabel("Frequency (MHz)")
    ax.set_ylabel("BP Coeff")

    title_str = ""
    if diff_thresh is not None:
        diff_str = "diff_thresh = %.2f" %(diff_thresh)
        title_str += diff_str
        if val_thresh is not None:
            title_str += ", "
        else: pass
    if val_thresh is not None:
        val_str = "val_thresh = %.2f" %(val_thresh)
        title_str += val_str
    else: pass

    if len(title_str):
        ax.set_title(title_str)
    else: pass

    ax.set_yscale('log')

    # If outfile, then save file, close window, and
    # turn interactive mode back on
    if outfile is not None:
        plt.savefig(outfile, dpi=150, bbox_inches='tight')
        plt.close()
    else:
        plt.show()

    return


def bp_filter(freqs, bp, diff_thresh=0.10, val_thresh=0,
              nchan_win=32, outfile=None):
    """
    Run the filter on a single bandpass file and find
    what channels need to be zapped

    diff_thresh = fractional diff threshold to mask chans
                  Mask if abs((bp-med)/med) > diff_thresh

    val_thresh  = min value threshold to mask chans
                  Mask if bp < val_thresh

    nchan_win = number of channels in moving window
    val_thresh  = min value threshold to mask chans
                  Mask if bp < val_thresh

    nchan_win = number of channels in moving window

    if outfile is specified, then save a plot showing
    the bandpass and masked channels
    """
    # Calculate running median and stdev
    mov_median, mov_std = moving_median(bp, nchan_win)

    # Calc fractional difference from median
    # Fix in case there are any zeros
    abs_med = np.abs(mov_median)
    if np.any(abs_med):
        eps = 1e-3 * np.min( abs_med[ abs_med > 0 ] )
    else:
        eps = 1e-3
    bp_diff = np.abs(bp - mov_median) / (abs_med + eps)

    # Find mask chans from diff
    diff_mask = np.where( bp_diff >= diff_thresh )[0]

    # Find mask chans from val
    bp_diff = np.abs(bp - mov_median) / (abs_med + eps)

    # Find mask chans from diff
    diff_mask = np.where( bp_diff >= diff_thresh )[0]

    # Find mask chans from val
    val_mask = np.where( bp < val_thresh )[0]

    # Get unique, sorted list of all bad chans
    mask_chans = np.unique( np.hstack( (diff_mask, val_mask) ) )

    # Get list of good chans (might need)
    all_chans = np.arange(0, len(freqs))
    good_chans = np.setdiff1d(all_chans, mask_chans)

    # Make a plot if outfile is specified
    if outfile is not None:
        plot_bp(freqs, bp, mask_chans, diff_thresh=diff_thresh,
                val_thresh=val_thresh, outfile=outfile)
    else:
        pass

    return mask_chans



############################
###   BANDPASS AND ZAP   ###
############################

def bp_zap(infile, png_dir, diff_thresh=0.25, nchan_win=5, 
           edgezap=0.15):
    """
    diff_thresh = fraction variation with running median 
                  bandpass to flag 
    nchan_win   = window size (in channels) for running median 
                  bandpass calculation 
    edge_zap    = zap out this fraction of the band at each edge
    """
    tstart = time.time()
    # get data
    tr = time.time()
    tt, freqs, dat = get_data(infile)
    dt_tr = time.time() - tr

    tc1 = time.time()
    if dat.dtype == 'float32':
        dat = dat.astype('float64')
    dt_tc1 = time.time() - tc1

    # calc bp 
    bp = np.mean(dat, axis=0)

    # find deviant channels
    tf = time.time()
    bchans = bp_filter(freqs, bp, diff_thresh=diff_thresh, 
                      nchan_win=nchan_win)
    dt_tf = time.time() - tf

    # calc edge chans
    nchans = len(bp)
    enum = int( edgezap * nchans )
    if enum <= 0:
        echans = np.array([])
    else:
        nn = np.arange(nchans)
        echans = np.hstack( (np.arange(0, enum+1), 
                            np.arange(nchans-enum, nchans)) )

    # Get unique set of mask chans
    mchans = np.unique( np.hstack( (bchans, echans) ) )

    outfile = infile.rsplit('.fil', 1)[0] + '_bp.png'
    basename = os.path.basename(outfile)
    outfile = png_dir + "/" + basename
    print(outfile)

    t_plt = time.time()
    plot_bp(freqs, bp, mchans, diff_thresh=diff_thresh, 
            val_thresh=None, outfile=outfile)
    dt_plt = time.time() - t_plt

    # Apply band pass and mean sub (ie, -1)
    t_ms = time.time()
    dd_corr = (dat / bp) - 1
    dt_ms = time.time() - t_ms

    # Calc channel norm 
    t_n = time.time()
    chan_norm = np.std(dd_corr, axis=0)
    t_n1 = time.time() 
    chan_norm[mchans] = 1
    t_n2 = time.time()
    dd_corr /= chan_norm
    t_n3 = time.time()
    dt_n = time.time() - t_n
    dt_n1 = t_n1 - t_n
    dt_n2 = t_n2 - t_n1
    dt_n3 = t_n3 - t_n2

    print(dd_corr.shape)

    # zero out zap chans
    t_z = time.time()
    dd_corr[:, mchans] = 0
    dt_z = time.time() - t_z

    tc2 = time.time()
    if dd_corr.dtype == 'float64':
        dd_corr = dd_corr.astype('float32')
    dt_tc2 = time.time() - tc2

    dt_tot = time.time() - tstart

    print("####################")
    print("##    BP ZAP      ##")
    print("####################")
    print(" read   = %.1f" %dt_tr)
    print(" dtype1 = %.1f" %dt_tc1)
    print(" dtype2 = %.1f" %dt_tc2)
    print(" plt    = %.1f" %dt_plt)
    print(" filter = %.1f" %dt_tf)
    print(" m sub  = %.1f" %dt_ms)
    print(" norm   = %.1f" %dt_n)
    print("   n1   = %.1f" %dt_n1)
    print("   n2   = %.1f" %dt_n2)
    print("   n3   = %.1f" %dt_n3)
    print(" total  = %.1f" %dt_tot)
    print("####################")

    return mchans, dd_corr


def combine_pol_dat(infile1, infile2, png_dir, diff_thresh=0.25, 
                    nchan_win=5, edgezap=0.15):
    """
    combine pol data
    """
    # Get mask chan and data for pol 1
    dtms = []
    for ii in range(5):
        tm = time.time()
        mm1, dd1 = bp_zap(infile1, png_dir, diff_thresh=diff_thresh, 
                          nchan_win=nchan_win, edgezap=edgezap) 
        dt_m = time.time() - tm
        dtms.append(dt_m)
    
    # Get full mask chans
    mm = np.unique( np.hstack( (mm1, mm1) ) )

    # Add together
    dd = (dd1 + dd1) / np.sqrt(2)

    # zap chans
    tz = time.time() 
    dd[:, mm] = 0
    dt_tz = time.time() - tz

    print("## COMBINE POLS ##")
    for ii, dt in enumerate(dtms):
        print(" mask%d - %.1f " %(ii, dt))
    print("##################")

    return dd


###########################
##   WRITING NEW FILE   ###
###########################

def write_fil(infile, outfile, dat, tel_id=None, 
              machine_id=None, src_name=None):
    """
    Write combined data (dat) to outfile.

    Use infile for header info (can be either pol) 
    """
    in_yr = your.Your(infile)

    if tel_id is not None:
        in_yr.telescope_id = tel_id

    if machine_id is not None:
        in_yr.machine_id = machine_id

    if src_name is not None:
        in_yr.your_header.source_name = src_name

    hdr_in = in_yr.your_header

    # Setup output filterbank file
    sig_obj = make_sigproc_object(
            rawdatafile   = outfile,
            source_name   = hdr_in.source_name,
            nchans        = hdr_in.nchans,
            foff          = hdr_in.foff,
            fch1          = hdr_in.fch1,
            tsamp         = hdr_in.tsamp,
            tstart        = hdr_in.tstart,
            src_raj       = in_yr.src_raj,
            src_dej       = in_yr.src_dej,
            machine_id    = in_yr.machine_id,
            nbeams        = in_yr.nbeams,
            ibeam         = in_yr.ibeam,
            nbits         = in_yr.nbits,
            nifs          = in_yr.nifs,
            barycentric   = in_yr.barycentric,
            pulsarcentric = in_yr.pulsarcentric ,
            telescope_id  = in_yr.telescope_id,
            data_type     = in_yr.data_type,
            az_start      = in_yr.az_start,
            za_start      = in_yr.za_start)

    # Write header
    sig_obj.write_header(outfile)

    # write data
    sig_obj.append_spectra(dat, outfile)

    return


def parse_input():
    """
    Use argparse to parse input
    """
    prog_desc = "BP correct, zap, and combine pols "+\
                "(can also change some header vals if you want)"
    parser = ArgumentParser(description=prog_desc)

    parser.add_argument('infile1', help='Pol 1 fil file')
    parser.add_argument('infile2', help='Pol 2 fil file')
    parser.add_argument('png_dir', help='Directory for png files')
    parser.add_argument('-tel', '--telescope_id',
                        help='Telescope ID Number',
                        required=False, type=int)
    parser.add_argument('-m', '--machine_id',
                        help='Machine ID Number',
                        required=False, type=int)
    parser.add_argument('-src', '--source_name',
                        help='Source Name',
                        required=False)
    parser.add_argument('-o', '--outbase',
                        help='Output file basename (no suffix)',
                        required=True)
    parser.add_argument('-ez', '--edgezap',
                        help='Fraction of band to zap at edges',
                        required=False, type=float, default=0.15)
    parser.add_argument('-dthresh', '--dthresh',
                        help='Fractional diff threshold for bp flagging',
                        required=False, type=float, default=0.25)
    parser.add_argument('-nwin', '--chanwin',
                        help='Window size for bandpass flagging (chans)',
                        required=False, type=int, default=5)

    args = parser.parse_args()

    return args



def combine_pol_fil(args):
    """
    Main function
    """
    infile1 = args.infile1
    infile2 = args.infile2

    tel_id = args.telescope_id
    machine_id = args.machine_id
    src_name = args.source_name

    edgezap = args.edgezap
    diff_thresh = args.dthresh 
    nchan_win = args.chanwin 

    outbase = args.outbase 
    outfile = outbase + '.fil'

    png_dir = args.png_dir

    dd = combine_pol_dat(infile1, infile2, png_dir, diff_thresh=diff_thresh,
                         nchan_win=nchan_win, edgezap=edgezap)

    tw = time.time()
    #write_fil(infile1, outfile, dd, tel_id=tel_id,
    #          machine_id=machine_id, src_name=src_name)
    dt_tw = time.time() - tw
    print("##  WRITE FIL  ##")
    print(" write - %.1f " %dt_tw)
    print("#################")

    return

################
###   MAIN   ###
################

debug = 0
if __name__ == "__main__":
    if not debug:
        args = parse_input()
        tstart = time.time() 
        combine_pol_fil(args)
        dt = time.time() - tstart
        print("Total Time: %.1f " %dt)
