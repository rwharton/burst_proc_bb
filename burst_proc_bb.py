import numpy as np
import subprocess
import time
import os
import sys
import shutil
import glob
import parse_arg as pa
import click

# set top of src directory (where this file is)
cur_dir = os.path.realpath(__file__)
srcdir  = cur_dir.rsplit('/', 1)[0]

def filter_cands(par):
    """
    Read in *singlepulse file 'spfile' and remove duplicate
    candidates, keeping the highest SNR version.
    """
    spfile = "%s/%s" %(par.outdir, par.sp_file)
    snr_min = par.sp_snr_min
    fac     = par.sp_match_fac

    if not os.path.exists(spfile):
        print("single pulse file: %s not found!" %spfile)
        return
    else: pass

    samples = []
    snrs    = []
    widths  = []
    lines   = []
    
    # read input file
    Ncand = 0
    with open(spfile, 'r') as fin:
        for line in fin:
            if line[0] in ["#", "\n"]:
                continue
            else: pass

            cols = line.split()
            
            if len(cols) not in [5, 6]:
                continue
            else: pass

            Ncand += 1
            snr   = float(cols[1])
            samp  = int(cols[3])
            width = int(cols[4])

            if snr < snr_min:
                continue
            else: pass

            samples.append( samp )
            widths.append( width )
            snrs.append( snr )
            lines.append(line)

    # remove dupes
    samples = np.array(samples)
    snrs = np.array(snrs)
    widths = np.array(widths)
    lines = np.array(lines)
    
    xx_nn = np.argsort(samples)
    samp_s  = samples[xx_nn]
    snr_s   = snrs[xx_nn]
    width_s = widths[xx_nn]

    uu_idx = []
    ii = 0
    while ii < len(xx_nn):
        samp_cur  = samp_s[ii]
        snr_cur   = snr_s[ii]
        width_cur = width_s[ii]

        samp_rest  = samp_s[ii+1:]
        snr_rest   = snr_s[ii+1:]
        width_rest = width_s[ii+1:]

        dn_abs = np.abs( samp_rest - samp_cur )
        dw = fac * 0.5 * ( width_cur + width_rest ) 

        xx_ii = np.where( dn_abs <= dw )[0]

        if len(xx_ii) > 0:
            jj = np.argmax( snr_rest[xx_ii] )
            if snr_rest[ xx_ii[jj] ] >= snr_cur:
                o_idx = ii + xx_ii[jj] + 1
            else:
                o_idx = ii
            ii += xx_ii[-1] + 2
        else:
            o_idx = ii
            ii += 1

        uu_idx.append(o_idx)

    idx = xx_nn[ uu_idx ]
    ulines = lines[ idx ]

    fbase, fsuf = spfile.rsplit('.', 1)
    outfile = "%s_unique.%s" %(fbase, fsuf)

    with open(outfile, 'w') as fout:
        for line in ulines:
            fout.write(line)
    
    print("\n\n=========== SP FILTER ==========")
    print("  input file:  %s" %spfile)
    print("  output file: %s" %outfile)
    print("  ")
    print(" Full cands      = %d" %Ncand)
    print(" SNR > %.1f      = %d" %(snr_min, len(lines)))
    print(" Filtered cands  = %d" %len(ulines))
    print("================================\n\n")

    out_sp_file = outfile.rsplit('/', 1)[-1]

    return out_sp_file
 

def vrad_2_cs(par):
    """
    Conversion from vrad to cs files
    """
    # Hold dict of params
    print("\n\n---------------------------------------------")
    print("Converting from vrad to cs")
    print("---------------------------------------------")
    params = {"inf": par.inf_file, "sp": par.sp_file, 
              "vrad_dir": par.vrad_dir, "nproc": par.nproc, 
              "vrad_base": par.vrad_base, "out": par.outdir, 
              "freq_dat": par.freq_dat, "freq_sp" : par.freq_sp, 
              "source": par.source, "tele": par.telescope, 
              "scan" : par.scan, "tskip": par.tskip, 
              "amount": par.data_amount, "dm": par.dm, 
              "src" : srcdir, "off_fn" : par.offset_file}

    # vrad2cs will make cs_dir and vdr_dir if they do not 
    # already exist.

    # NOTE: 'inf_file' and 'sp_file' need to be in the 
    #        output directory
    
    cmd = "python3 -u %(src)s/vrad2cs/vrad2cs.py" \
           " --inf_file %(out)s/%(inf)s" \
           " --off_file %(off_fn)s" \
           " --sp_file %(out)s/%(sp)s --vrad_dir %(vrad_dir)s" \
           " --vrad_base %(vrad_base)s --out_dir %(out)s" \
           " --freq_band %(freq_dat)s --sp_freq_band %(freq_sp)s" \
           " --nproc %(nproc)d --dm %(dm)f --source %(source)s" \
           " --scan %(scan)d --tskip %(tskip)d" \
           " --telescope %(tele)s --data_amount %(amount)i" % params
    
    print(cmd)
    subprocess.run(cmd, shell=True)

    return


def cs_2_fil(par):
    """
    Converts cs to filterbank for all the channel numbers given
    """
    print("\n\n---------------------------------------------")
    print("Converting from cs to fil")
    print("---------------------------------------------")
    
    params = {"out": par.outdir, "base": par.cs_base, 
              "cs_dir": par.cs_dir, "dada_dir": par.dada_dir, 
              "fil_dir": par.fil_dir, "dm": par.dm, 
              "src" : srcdir, "nproc": par.nproc}

    # cs2fil will make fil_dir and dada_dir if they don't exist

    for num in par.nchans:    
        params["num"] = num
        cmd = "python3 -u %(src)s/cs2fil/cs2fil.py" \
              " --basename %(base)s --cs_dir %(out)s/%(cs_dir)s"\
              " --dada_dir %(out)s/%(dada_dir)s" \
              " --fil_dir %(out)s/%(fil_dir)s --dm %(dm)f" \
              " --nchan %(num)d --nproc %(nproc)d" % params
        
        print(cmd)
        subprocess.run(cmd, shell=True)

    return


def combine_pol(par):
    """
    Combines lcp and rcp files
    """
    print("\n\n---------------------------------------------")
    print("Combining polarizations")
    print("---------------------------------------------")
    # Go through each channel number
    params = {"out": par.outdir, "lcp_base": par.lcp_base, 
              "rcp_base": par.rcp_base, 
              "comb_base": par.comb_base, "fil_dir": par.fil_dir, 
              "png_dir": par.png_dir, "ez": par.ez, "dthresh": par.dthresh,
              "nwin": par.nwin, "src" : srcdir, "nproc" : par.nproc}
    
    cmd = "python3 -u %(src)s/combine_pol/combine_pol.py" \
          " --lcp_base %(lcp_base)s --rcp_base %(rcp_base)s" \
          " --comb_base %(comb_base)s --fil_dir %(out)s/%(fil_dir)s" \
          " --png_dir %(out)s/%(png_dir)s --edgezap %(ez)f" \
          " --dthresh %(dthresh)f --chanwin %(nwin)i" \
          " --nproc %(nproc)d" % params
        
    print(cmd)
    subprocess.run(cmd, shell=True)

    return


def plot_fil(par):
    """
    Plots filterbank files
    """
    print("\n\n---------------------------------------------")
    print("Converting from fil to plot")
    print("---------------------------------------------")

    # Check if png dir exists, make it if not
    full_png_dir = "%s/%s" %(par.outdir, par.png_dir)
    if not os.path.exists(full_png_dir):
        os.makedirs(full_png_dir)

    params = {"out": par.outdir, "offset": par.offset_file, 
              "fil_base": par.fil_base, "sp_dt" : par.sp_dt_us, 
              "fil_dir": par.fil_dir, "png_dir": par.png_dir,
              "dm": par.dm, "plt_dt": par.plt_dt, 
              "plt_chans" : par.plt_chans, "tdur": par.tdur, 
              "comb_base":par.comb_base, "src" : srcdir, 
              "nproc":par.nproc}

    cmd = "python3 -u %(src)s/plotfil/plotfil.py" \
            " --off_file %(out)s/%(offset)s --sp_dt %(sp_dt).4f" \
            " --fil_dir %(out)s/%(fil_dir)s --fil_base %(fil_base)s" \
            " --comb_base %(comb_base)s --png_dir %(out)s/%(png_dir)s" \
            " --dm %(dm)f --plt_dt %(plt_dt)f --tdur %(tdur)f" \
            " --plt_chans %(plt_chans)d --nproc %(nproc)d" % params
    
    print(cmd)
    subprocess.run(cmd, shell=True)
    
    return


def dm_opt(par):
    """
    Optimizes dm
    """
    print("\n\n---------------------------------------------")
    print("Finding optimal dm value")
    print("---------------------------------------------")
    params = {"out_dir": par.outdir, "dm_dir": par.dm_dir, 
              "fil_dir": par.fil_dir, "offset": par.offset_file, 
              "dm_lo": par.dm_lo, "dm_hi": par.dm_hi, 
              "dm_step": par.dm_step, "src" : srcdir}
    
    cmd = "python3 -u %(src)s/dm_opt/dm_opt.py --out_dir %(out_dir)s"\
           " --dm_dir %(dm_dir)s --fil_dir %(fil_dir)s"\
           " --offsets %(offset)s --dm_lo %(dm_lo)f" \
           " --dm_hi %(dm_hi)f --dm_step %(dm_step)f" % params
   
    print(cmd)
    subprocess.run(cmd, shell=True)

    return


def plot_scint(par):
    """
    Plots scintilation bandwidth
    """
    print("\n\n---------------------------------------------")
    print("Plotting Scintillation Bandwidth")
    print("---------------------------------------------")

    # Check if png dir exists, make it if not
    full_png_dir = "%s/%s" %(par.outdir, par.png_dir)
    if not os.path.exists(full_png_dir):
        os.makedirs(full_png_dir)

    params = {"out": par.outdir, "offset": par.offset_file, 
              "fil_dir": par.fil_dir, "png_dir": par.png_dir, 
              "dm": par.dm, "tdur": par.tdur, "nproc": par.nproc,
              "sp_dt" : par.sp_dt_us, "comb_base":par.comb_base, 
              "src" : srcdir}
    
    cmd = "python3 -u %(src)s/plot_scint/plot_scint.py" \
            " --off_file %(out)s/%(offset)s"\
            " --fil_dir %(out)s/%(fil_dir)s --comb_base %(comb_base)s "\
            " --png_dir %(out)s/%(png_dir)s --sp_dt %(sp_dt).4f" \
            " --dm %(dm)f --tavg 1 --tdur %(tdur)f" \
            " --nproc %(nproc)d" % params
    
    print(cmd)
    subprocess.run(cmd, shell=True)
    
    return


def get_pars(parfile):
    """
    Read in the parameter file, parse it, and 
    return a class called 'par'

    This is a way to simulate what we did before 
    of simply importing the parameter file as a 
    module.
    """
    pdict = pa.get_par_dict(parfile)
    par = pa.Parameters(pdict)
    return par


def print_times(times):
    """
    Times is a list of [lable, time] tuples 

    print summary
    """
    print("\n\n")
    print("################################")
    print("###      TIME SUMMARY        ###")
    print("################################")
    for tt in times:
        print("{:<16s}:  {:>12.2f}s".format(tt[0], tt[1]))
    print("################################")
    print("\n\n")
    return
    


@click.command()
@click.argument("parfile", type=str)
def run_pipeline(parfile):
    """
    Run the pipeline
    """
    # Get pars 
    par = get_pars(parfile)
    
    print('vrad', par.vrad_to_cs)
    print('cs', par.cs_to_fil)
    print('plotfil', par.plot_fil)
    print('plot_scint', par.plot_scint)
    print('combine', par.combine_pol)

    times = []

    if par.sp_filter:
        st = time.time()
        filter_file = filter_cands(par)
        par.sp_file = filter_file
        times.append(["filter dupes", round(time.time()-st, 2)])

    if par.vrad_to_cs:
        st = time.time()
        vrad_2_cs(par)
        times.append(["vrad->cs", round(time.time()-st, 2)])
    
    if par.cs_to_fil:
        st = time.time()
        cs_2_fil(par)
        times.append(["cs->fil", round(time.time()-st, 2)])
    
    if par.combine_pol:
        st = time.time()
        combine_pol(par)
        times.append(["Combine pols", round(time.time()-st, 2)])
    
    if par.plot_fil:
        st = time.time()
        plot_fil(par)
        times.append(["fil plot", round(time.time()-st, 2)])
    
    if par.plot_scint:
        st = time.time()
        plot_scint(par)
        times.append(["scint plot", round(time.time()-st, 2)])
    
    if par.dm_opt:
        st = time.time()
        dm_opt(par)
        times.append(["dm optimize", round(time.time()-st, 2)])
   
    total = sum(item[1] for item in times)
    times.append(["total", total])
    print_times(times)
    
    # clean up dada and vdr files
    if par.cleanup:
        dada_path = "%s/%s" % (par.outdir, par.dada_dir)
        if os.path.exists(dada_path):
            files = glob.glob(dada_path + '/*')
            for dfile in files:
                os.remove(dfile)
        vdr_path = "%s/%s" % (par.outdir, par.vdr_dir)
        if os.path.exists(vdr_path):
            files = glob.glob(vdr_path + "/*")
            for vfile in files:
                os.remove(vfile)
    return


debug = 0

if __name__ == "__main__":
    if not debug: 
        run_pipeline()
