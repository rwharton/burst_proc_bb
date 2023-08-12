# burst_proc_bb
Software to extract baseband data around candidate FRBs, 
create filterbanks, and make various diagnostic plots.

Designed to work (more or less) exclusively with DSN data 
products and assuming the JPL singularity environment.

The first version of the pipeline was written by Evan Zhang (Caltech) 
using some scripts written by me and Aaron Pearlman (McGill).  Final 
version (and all bugs) by me.

## Processing Options

Available processing options are:

  * Filter single pulse candidates 
  * Convert vrad file to cs
  * Produce filterbank(s) from cs
  * Bandpass, RFI flag, and combine pols
  * Plot burst dynamic spectrum
  * Estimate scintillation bandwidth and plot ACF
 
## Running 

Running is (quite possibly) very easy.  You just run 
the script and pass all the parameters through a parameter 
file:

    python burst_proc_bb.py params.txt

The parameter file is a text file and is parsed by the 
script.  It is *not* a python module that is imported. 
This is only important in so far as you cannot use other 
variable names or math / string operations to construct new 
ones.  Every variable and value are read simply by splitting 
the line at the equals sign and taking a key and value. 

## Parameter File

Here we will go through an example parameter file, which can 
be found as `params/params_23m041_sband.txt`.  We will first 
provide a section of the parameter file and then explain it.

    ##########################
    ###  INPUT VRAD FILES  ###
    ##########################
    
    # Directory containing the vrad files
    vrad_dir = "/raw/23m041"
    
    # vrad base
    vrad_base = "23-041"
    
    # Names of info / sp time files
    #  -->  Place these in output directory
    inf_file = "scan.table.23m041"
    sp_file  = "slcp.singlepulse"
    
    # Filter sp cand file to remove duplicates?
    # Can also set SNR threshold and sp_match_fac,
    # which grows the min separation by which cands
    # are said to be matches (1.0 means just the widths)
    sp_snr_min   = 10
    sp_match_fac = 1.25
    
    # Time resolution (in us) of singlepulse file
    sp_dt_us = 10.24
   
 
In the first part, you need to specify where the raw data 
(`.vrad`) files are using `vrad_dir`.  The `vrad_base` parameter 
is the base name of the `.vrad` files.  Using this basename and 
the frequency band (provided later), the script will glob on 
the base name in the vrad directory.  

You also need to give the names of files providing meta data info 
and the arrival time of the burst candidates.  The `inf_file` is 
either a PRESTO `.inf` file or a scan table file.  The `sp_file` is 
the PRESTO singlepulse file.  These files both need to be in the output 
directory (specified later).  The input single pulse file can 
be filtered with a minimum SNR cutoff (`sp_snr_min`) and by 
removing duplicates found within `sp_match_fac` times the width 
of a candidate.  Finally, the time resolution used for the single 
pulse search is needed (`sp_dt_us`) to convert pulse widths in bins 
to seconds.  The time resolution should be specified here in microseconds. 
    
    ############################
    ###  OUTPUT DIRECTORIES  ###
    ############################
    
    # Output directory
    outdir = "/data/frb_pipe/23m041/sband"
    
    # Output subdirectories (shouldn't need to change)
    cs_dir = "cs"
    dada_dir = "dada"
    fil_dir = "fil"
    png_dir = "png"
    vdr_dir = "vdr" 

These parameters specify the output directories.  You should only 
ever change `outdir`, which specifies the top of the output directory. 
This is where the scan table and single pulse files should go.  The 
subdirectories will organize the various data products.  Don't change these.

    ######################
    ###  STEPS TO RUN  ###
    ######################
    
    # Use (0 or False) and (1 or True) to
    # set which steps to run
    
    sp_filter   = 1
    vrad_to_cs  = 1
    cs_to_fil   = 1
    combine_pol = 1
    plot_fil    = 1
    plot_scint  = 1
    
    # max number of processes to run at
    # once with parallel processing
    nproc = 4

In this section we set which steps we want to run.  You can either 
use (0 = dont run, 1 = run) or (False = don't run, True = run) to 
select which you want to run.  Most steps will require some previous 
data product to be available.  It will check that they exist first, 
and skip otherwise.  The `nproc` parameter sets the number of processes 
to run in parallel. 

    ######################
    ###   VRAD -> CS   ###
    ######################
    
    # Frequency configuration of data being processed (freq_dat)
    # and of the single pulse file (freq_sp).  The names for these
    # configurations are defined in vrad2cs/freq_setups.py
    freq_dat = 'sband-1'
    freq_sp  = 'sband-1'
    
    source = "J1810-197"
    telescope = "robledo"
    data_amount = 3 # in seconds

These parameters deal with the conversion from the telescope 
data format (vdr) to a format readable by SIGPROC and other 
pulsar software (cs).  The `source` and `telescope` specify 
data about the observation and `data_amount` specifies how much 
data should be read around the burst (needs to be integer second).

The first two parameters (`freq_dat` and `freq_sp`) give the frequency 
setup for the data we are processing (`freq_dat`) and the data that 
was used for the initial single pulse search (`freq_sp`).  These may 
be different when, for example, you want to extract X-band data around 
S-band bursts.  The frequency setups are named and specified in 
`vrad2cs/freq_setups.py` under `get_freq_info()`:

    if fkey == "sband-1":
        fsetup.name      = fkey
        fsetup.band      = "s"
        fsetup.pol_names = ["SLCP", "SRCP"]
        fsetup.sub_freqs = [2250]
        fsetup.sub_bw    = 100

    elif fkey == "xband-1":
        fsetup.name      = fkey
        fsetup.band      = "x"
        fsetup.pol_names = ["XLCP", "XRCP"]
        fsetup.sub_freqs = [8224, 8256, 8288, 8320]
        fsetup.sub_bw    = 32

If you want to add any other frequency setups, you can do that there.

    #####################
    ###   CS -> FIL   ###
    #####################
    
    cs_base = "s"
    dm = 178.850
    nchans = [8000] 

These parameters are used during the conversion from complex sampled 
baseband data (cs) to filterbank.  The `cs_base` is used to glob 
on the needed cs files, which will generally start with (slcp/srcp 
or xlcp/xrcp).  In this case where we have s-band data, globbing on 
"s" is sufficient. 

The `dm` here is the dispersion measure of the burst, which will be 
used for coherent de-dispersion when making the filterbank and used 
for de-dispersing in plots later.  The `nchans` parameter is the 
number of channels to make in the filterbank.  This has to be a list, 
but you don't have to give more than one value if you don't want.  If 
you do give multiple values, then this will produce multiple filterbanks.


    ##############################
    ###   ZAP + COMBINE POLS   ###
    ##############################
    
    lcp_base = "slcp"
    rcp_base = "srcp"
    comb_base = "sbnd"
    ez = 0.10  # edge zap
    dthresh = 10.0 # Fractional diff threshold
    nwin = 5 # window size

These parameters are for the rfi flaggin, bandpass correction, 
and polarization summing steps.  The base names give the first 
parts of the filterbank files.  Generally these are just 
"slcp/srcp" or "xlcp/xrcp".  The remaining parameters are 
the edge zapping (`ez`) and the threshold to zap bad channels 
before combining the polarizations.

    ######################
    ###   PLOT BURST   ###
    ######################
    
    offset_file = "offset_times_sband.txt"
    fil_base  = "s"
    plt_dt    = 0.1  # plot at this time res in ms
    plt_chans = 100  # plot at most this many chans
    tdur = 0.2       # grab this much data for plot

These parameters are for the plotting of burst candidates.
The offset file contains the burst time relative to the 
start of the (usually) 3-sec data chunk. The `fil_base` 
is the base name of the fil files that will be globbed for. 
By using just "s" or "x" you can get both pols and the 
combined data.  If you just want combined then you could 
do "sbnd" or "xbnd".  The `plot_dt` and `plt_chans` give 
the time resolution and channel number for the plots.  This 
is useful because the intrinsic res will likely be much 
finer than you'd like for a plot.  The `tdur` parameter 
sets how much data to read for plotting, you shouldn't really 
need to change that.

    ###########################
    ###   DM OPTIMIZATION   ###
    ###########################
    
    # Search over trial DMs to find max snr
    # probably not that useful
    
    dm_opt = False   # run optimization?
    dm_lo = 0
    dm_hi = 400
    dm_step = 20
    dm_dir = "dm"

This will be for DM optimization, but is not currently implemented.
Don't need to set these.

    ####################
    ###   CLEAN UP   ###
    ####################
    
    # Cleanup (remove dada and vdr files)
    cleanup = True

Remove intermediary data products.
            
