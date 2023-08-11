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
    
    
    
