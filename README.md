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
 

