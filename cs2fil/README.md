# cs2fil
Preliminary processing of DSN baseband data

    Usage: cs2fil.py [OPTIONS]
    
      Convert multiple chunks of complex sampled voltage data  to coherently de-
      dispersed channelized filterbanks.
    
      Input cs files assumed to be in the form
    
             {basename}-XXXX-YY.cs in directory
    
      where {basename} is the base name associated with all  files, XXXX is a
      number 0000-9999 that gives the data  chunk order, and YY is 01-07 and
      represents the channel.
    
      The channels of each time chunk will be combined into a  DADA file called
    
             {basename}-XXXX.dada
    
      Then run digifil to produce a coherently de-dispersed  filterbank at DM "dm"
      with a total of nchan channels. Use nthread threads for digifil and remove
      incohrent  dispersive delay if inc_ddm=True
    
      The filterbank file will be called
    
            {basename}-XXXX.fil
    
    Options:
      --basename TEXT    Name of cs files: {basename}*.cs
      --cs_dir TEXT      Directory containing cs files
      --dada_dir TEXT    Directory for DADA files
      --fil_dir TEXT     Directory for *.fil files
      --dm FLOAT         DM for intra-channel coherent dedispersion
      --nchan INTEGER    Total number of output channels in filterbank
      --mem_lim FLOAT    Max memory to use during DADA conversion (GB)
      --nthread INTEGER  Number of threads for digifil processing
      --help             Show this message and exit.
