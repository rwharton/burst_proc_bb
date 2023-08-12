import numpy as np
import your
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


# Lorimer & Kramer (2006)
kdm = 4148.808 # MHz^2 / (pc cm^-3)

def dm_delay(f1, f2, DM, kdm=kdm):
    """
    Returns DM delay in sec
    """
    return kdm * DM * (1.0 / (f1 * f1) - 1 / (f2 * f2))


def deltaT_old(ichan, dt, df, f0, kdm=kdm):
    """
    Return DM=1pc/cc delay (in samples) for channel number ichan

    dt = time resolution in sec
    df = channel width in MHz
    f0 = first channel frequency in MHz
    ichan = channel number (0 index)
    """
    return (kdm / dt) * ((f0 + ichan * df)**-2.0 - f0**-2)


def dmdt_old(DM, ichan, dt, df, f0, kdm=kdm):
    """
    Return DM delay for DM=DM in *integer* samples
    """
    return int(np.round( DM * deltaT(ichan, dt, df, f0, kdm=kdm)))


def dmdt_float_old(DM, ichan, dt, df, f0, kdm=kdm):
    """
    Return DM delay in *float* samples
    """
    return DM * deltaT(ichan, dt, df, f0, kdm=kdm)


def deltaT(freqs, f0, dt, kdm=kdm):
    """
    Return array of single unit DM delays in samples for MHz freqs
    """
    return (kdm / dt) * (freqs**-2.0 - f0**-2)


def dmdt(DM, freqs, f0, dt, kdm=kdm):
    """
    Return DM delays (in samples) for DM=DM 
    """
    return DM * deltaT(freqs, f0, dt, kdm=kdm)


def dedisperse_one(dspec, dm, dt, df, f0, kdm=kdm):
    """
     
    """
    nchans = dspec.shape[0]
    nt = dspec.shape[1]
    dsamps = np.array([ dmdt_old(dm, ichan, dt, df, f0, kdm=kdm) for ichan in range(nchans) ])
    dsamps -= np.min(dsamps)
    tpad = np.max(dsamps)
    outarr = np.zeros( nt + tpad )
    for ii in range(nchans):
        osl = slice(tpad - dsamps[ii], nt + tpad - dsamps[ii])
        outarr[osl] += dspec[ii]
    return outarr[tpad:nt + tpad] / float(nchans)


def dedisperse_dspec_old(dspec, dm, dt, df, f0, kdm=kdm):
    nchans = dspec.shape[0]
    nt = dspec.shape[1]
    dsamps = np.array([ dmdt_old(dm, ichan, dt, df, f0, kdm=kdm) for ichan in range(nchans) ])
    dsamps -= np.min(dsamps)
    tpad = np.max(dsamps)
    outarr = np.zeros( (nchans, nt + tpad) )
    for ii in range(nchans):
        osl = slice(tpad - dsamps[ii], nt + tpad - dsamps[ii])
        outarr[ii, osl] = dspec[ii]
    return outarr[:, tpad:nt + tpad]


def dedisperse_dspec(dspec, dm, freqs, f0, dt, kdm=kdm, reverse=False):
    nchans = dspec.shape[0]
    nt = dspec.shape[1]
    dsamps = dmdt(dm, freqs, f0, dt, kdm=kdm)
    #dsamps -= np.min(dsamps)
    if reverse:
        sgn = -1.0
    else:
        sgn = +1.0 
   
    nt_out = nt
    if nt_out % 2:
        nt_out -= 1
    else: pass

    dout = np.zeros( (nchans, nt_out) )
    
    for ii, dd in enumerate(dspec):
        ak = np.fft.rfft(dd)
        bfreq = np.arange(len(ak)) / (1.0 * len(dd))
        shift = np.exp(sgn * 1.0j * 2 * np.pi * bfreq * dsamps[ii])
        dd_shift = np.fft.irfft( ak * shift )
        #print(len(dd), len(ak), len(dd_shift))
        dout[ii] = dd_shift[:]

    return dout


def dspec_avg_chan(dspec, freqs, avg_chan=1):
    Nchan = dspec.shape[1]
    n = int(Nchan / avg_chan)

    freq_out = np.zeros(n)
    dd_out = np.zeros( (dspec.shape[0], n) )

    for ii in range(n):
        sl = slice(ii * avg_chan, (ii+1) * avg_chan)
        freq_out[ii] = np.mean(freqs[sl])
        dd_out[:, ii] = np.mean(dspec[:, sl], axis=1)

    return freq_out, dd_out


def dspec_avg_time(dspec, tt, avg_tsamp=1):
    Nt = dspec.shape[0]
    n = int(Nt / avg_tsamp)

    tt_out = np.zeros(n)
    dd_out = np.zeros( (n, dspec.shape[1]) )

    for ii in range(n):
        sl = slice(ii * avg_tsamp, (ii+1) * avg_tsamp)
        tt_out[ii] = np.mean(tt[sl])
        dd_out[ii, :] = np.mean(dspec[sl,:], axis=0)

    return tt_out, dd_out


def get_chan_info(data_file):
    """
    Get obs info
    """
    yr = your.Your(data_file)
    foff = yr.your_header.foff
    fch1 = yr.your_header.fch1
    dt   = yr.your_header.tsamp
    nchans = yr.your_header.nchans
    return nchans, fch1, foff, dt


def get_snippet_data(filfile, dm, favg=1, tavg=1, bpass=True,
                     tp=None, tdur=None):
    """
    Use your to read data and get metadata
    """
    # Get info
    nchans, fch1, foff, dt = get_chan_info(filfile)
    yr = your.Your(filfile)
    nsamps = yr.your_header.nspectra
    freqs = np.arange(nchans) * foff + fch1
    tt = np.arange(nsamps) * dt
    tt -= np.mean(tt)
    
    # Get data
    # if pulse start + dur not stated, read full file
    if (tp is None) or (tdur is None):
        dat = yr.get_data(0, yr.your_header.nspectra)
    else:
        nsamps = int( tdur / dt )
        npeak  = int( tp / dt )
        nstart = npeak - nsamps // 2
        dat = yr.get_data(nstart, nsamps)
        tt = np.arange(nsamps) * dt
        tt -= np.mean(tt)

    # Dedisperse
    dout = dedisperse_dspec(dat.T, dm, freqs, freqs[0], dt) 
    dout = dout.T
    if bpass:
        Nthird = int( nsamps / 3 )
        bp = np.mean(dout[:Nthird], axis=0)
        xx = np.where(bp==0)[0]
        if len(xx):
            bp[xx] = 1
        dout = dout/bp - 1
    else: pass

    # Average in time (if desired)
    if tavg > 1:
        tt_out, dout = dspec_avg_time(dout, tt, avg_tsamp=tavg)
    else:
        tt_out = tt

    # Average in freq (if desired)
    if favg > 1:
        ff_out, dout = dspec_avg_chan(dout, freqs, avg_chan=favg)
    else:
        ff_out = freqs

    return tt_out, ff_out, dout, dt


def make_plot(filfile, dm, pulse_num, snr, 
              plt_chans=100, plt_dt=1.0, spec_sig=5,
              outbins=128, tp=None, tdur=None, 
              outfile=None, cnum=None, ctime=None, 
              cdm=None, cwidth=None, pulse_width=None, 
              tcand=None, bpass=True):
    """
    make 3 panel plot
    """
    if outfile is not None:
        plt.ioff()
    else: pass 

    nchan0, _, _, dt0 = get_chan_info(filfile)

    tavg = max( 1 , int(plt_dt*1e-3 / dt0) ) 
    favg = max( 1 , int(nchan0 / plt_chans) )

    tt, freqs, dat, dt = get_snippet_data(filfile, dm, 
                             favg=favg, tavg=tavg, bpass=bpass,
                             tp=tp, tdur=tdur)
    tt *= 1e3
    #tt0, _, dat0 = get_snippet_data(filfile, 0, 
    #                         favg=favg, tavg=tavg, bpass=True)
    #tt0 *= 1e3

    # Set tlim if desired
    if outbins > 0:
        Nt = len(tt)
        xt_lo = max( 0, Nt//2 - outbins//2)
        xt_hi = min( Nt-1, Nt//2 + outbins//2)
        tlim = (tt[xt_lo], tt[xt_hi])
    else:
        tlim = (tt[0], tt[-1])

    fig = plt.figure(constrained_layout=True)
    gs = GridSpec(3, 3, figure=fig)

    ax_t  = fig.add_subplot(gs[0, 0:2])
    ax_ds = fig.add_subplot(gs[1:, 0:2])
    ax_f  = fig.add_subplot(gs[1:, 2])
    ax_txt = fig.add_subplot(gs[0, 2])


    # Dynamic Spectrum
    ext = [tt[0], tt[-1], freqs[0], freqs[-1]]
    d_sig = np.std(dat)
    d_med = np.median(dat)
    vmax = d_med + 4 * d_sig
    vmin = d_med - 3 * d_sig
    ax_ds.imshow(dat.T, aspect='auto', interpolation='nearest',
                 origin='lower', extent=ext, vmin=vmin, vmax=vmax)
    ax_ds.set_xlabel("Time (ms)", fontsize=14)
    ax_ds.set_ylabel("Freq (MHz)", fontsize=14)
    ax_ds.set_xlim(tlim)

    # Time series
    ts = np.mean(dat, axis=1)
    Nthird = len(ts) // 3
    avg = np.mean(ts[:Nthird])
    sig = np.std(ts[:Nthird])
    ts = (ts-avg)/sig
    
    #ts0 = np.mean(dat0, axis=1)
    #ts0 = (ts0-avg)/sig
    #ax_t.plot(tt0, ts0, c='0.5', zorder=-1)
    ax_t.plot(tt, ts)
    ax_t.tick_params(axis='x', labelbottom=False)
    tylim = ax_t.get_ylim()
    #ax_t.set_xlim(tt[0], tt[-1])
    ax_t.set_xlim(tlim)

    # Spectrum
    xpk = np.argmax(ts)
    #xx_below = np.where( ts <= 0.1*np.max(ts) )[0]
    #xx_lo = np.max(xx_below[xx_below <= xpk ])
    #xx_hi = np.min(xx_below[xx_below >= xpk ])

    xx_lo = int( len(ts) // 2 ) - 2
    xx_hi = int( len(ts) // 2 ) + 2

    if pulse_width:
        sample_width = max(pulse_width/dt/tavg, 2)
        xx_lo = int( len(ts) // 2 ) - int(sample_width//2)
        xx_hi = int( len(ts) // 2 ) + int(sample_width//2)
    
    off_spec = np.mean(dat[:Nthird], axis=0)
    off_sig = np.std(off_spec)
    off_sig = off_sig * np.sqrt(Nthird/(xx_hi-xx_lo))

    spec = np.mean(dat[xx_lo:xx_hi], axis=0) / off_sig
    ax_f.plot(spec, freqs)
    ax_f.set_ylim(freqs[0], freqs[-1])
    ax_f.tick_params(axis='y', labelleft=False)

    # Shade region used for spec
    tlo = tt[xx_lo]
    thi = tt[xx_hi]
    ax_t.fill_betweenx([-10, 100], tlo, thi, color='r', alpha=0.1)
    ax_t.set_ylim(tylim)

    # Add cand info subplot
    ax_txt.axis('off')
    outstr = ""
    if cnum is not None:
        outstr += "Cand: %d\n" %cnum 
    if ctime is not None:
        outstr += "Time: %.3f s\n" %ctime
    if cdm is not None:
        outstr += "DM: %.2f\n" %cdm
    if cwidth is not None:
        outstr += "Width: %d bins\n" %cwidth
    ax_txt.text(0.00, 0.8, outstr, fontsize=12, ha='left', va='top', 
                transform=ax_txt.transAxes)
    
    # Add pulse num and snr text
    if pulse_width and tcand:
        txt_str = "Pulse: %d\n" %pulse_num +\
                  "SNR: %.1f\n" %snr +\
                  "Width: %.3f ms\n" %(pulse_width * 1e3) +\
                  "DM: %.1f pc/cc\n" %dm +\
                  "tcand: %.4f\n" %tcand
        fig.text(0.75, 0.95, txt_str, fontsize=12, ha='left', va='top')
    else:
        txt_str = "Pulse: %d\n" %pulse_num +\
                  "SNR: %.1f\n" %snr +\
                  "DM: %.1f pc/cc\n" %dm
        fig.text(0.75, 0.95, txt_str, fontsize=12, ha='left', va='top')

    if outfile is not None:
        plt.savefig(outfile, dpi=100, bbox_inches='tight')
        plt.close()
        plt.ion()
    else: 
        plt.show()

    return freqs, spec


