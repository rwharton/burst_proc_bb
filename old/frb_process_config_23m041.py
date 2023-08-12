"""
Parameter file for frb_process.py
"""


# Input directory
indir = "/src/in/"

# Output directory
outdir = "/out/23m041/"

# Inputs for vrad2cs conversion
inf_file = "scan.table.23m041"
sp_file = "one.singlepulse"
vrad_dir = "vrad/"
vrad_base = "23-041" # first 5 digits identifying telescope reading

freq_band = "x" # set to x or s 
source = "j1810" 
telescope = "robledo" 
data_amount = 3 # in seconds

# Inputs for cs2fil conversion
cs_base = "x" # set to x or s if you want same files to be convereted
cs_dir = "cs/"
dada_dir = "dada/"
fil_dir = "fil/"
dm = 178.85
nchans = [100] # list

# Inputs for combine_pol
lcp_base = "xlcp"
rcp_base = "xrcp"
comb_base = "xbnd"
ez = 0.15 # edge zap
dthresh = 0.25 # Fractional diff threshold
nwin = 5 # window size

# inputs for plotfil
offset_file = "offset_times.txt"
fil_base = "x" # set to x or s if you want same files to be converted
png_dir = "png/"
npy_dir = "npy/"
tavg = 100
tdur = 0.2

# scintilation inputs
npy_base = "x"

# Conversions to run
vrad_to_cs = True
cs_to_fil = True
plot_fil = True
plot_scint = True
combine_pol = True

# Optimizations to run
dm_opt = True

# DM Optimizition
dm_lo = 0
dm_hi = 400
dm_step = 100
dm_dir = "dm/"

# Cleanup (remove dada and vdr files)
cleanup = True
vdr_dir = "vdr/"