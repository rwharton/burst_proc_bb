import numpy as np
import your

####################
##  GETTING DATA  ##
####################

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


def get_spec(filfile, dm, favg=1, tavg=1, tp=None, tdur=None, 
             pulse_width=None, bpass=True):
    """
    Get spectrum from pulse filterbank
    """
    tt, freqs, dat, dt = get_snippet_data(filfile, dm, 
                             favg=favg, tavg=tavg, bpass=bpass,
                             tp=tp, tdur=tdur)
    tt *= 1e3
    
    # Time series
    ts = np.mean(dat, axis=1)
    Nthird = len(ts) // 3
    avg = np.mean(ts[:Nthird])
    sig = np.std(ts[:Nthird])
    ts = (ts-avg)/sig
    
    # Spectrum
    xpk = np.argmax(ts)
    xx_lo = int( len(ts) // 2 ) - 2
    xx_hi = int( len(ts) // 2 ) + 2

    if pulse_width:
        sample_width = pulse_width/dt/tavg
        xx_lo = int( len(ts) // 2 ) - int(sample_width//2)
        xx_hi = int( len(ts) // 2 ) + int(sample_width//2)
    
    off_spec = np.mean(dat[:Nthird], axis=0)
    off_sig = np.std(off_spec)
    off_sig = off_sig * np.sqrt(Nthird/(xx_hi-xx_lo))

    spec = np.mean(dat[xx_lo:xx_hi], axis=0) / off_sig

    return freqs, spec


