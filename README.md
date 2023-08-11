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
be found as `params/params_23m041_sband.txt`.

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
    
    
    
    
