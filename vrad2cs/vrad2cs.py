import numpy as np
import os 
import sys
import glob
import time
import subprocess
import click
from datetime import datetime, timedelta
from astropy.time import Time, TimeDelta
import sigproc as fb
import multiprocessing as mp
import re
import freq_setups as fs

cur_dir = os.path.realpath(__file__)
srcdir  = cur_dir.rsplit('/', 1)[0]

def extract_start_time(inf_file):
    """
    Given inf_file/scan_table, find MJD time and return it 
    in YY/Day_HH:MM:SS format
    """
    if "scan.table" in inf_file:
        # Define the output format
        output_format = "%y/%j_%H:%M:%S"
        # Open the scan table file
        with open(inf_file, 'r') as file:

            line = file.readline().rstrip()

            # Split the line into parts
            parts = line.split()

            # Extract the date and time parts
            date = parts[3] # could be different
            time = parts[4] + parts[5] + parts[6] # could be different

            date_time_str = date + time.replace(' ', '')
            print(date_time_str)
            date_time_obj = datetime.strptime(date_time_str, '%y%m%d%H%M%S')
            output = date_time_obj.strftime(output_format)
            print(Time(date_time_obj).yday)
            return Time(date_time_obj), output
    else:
        epoch_key = "Epoch of observation (MJD)"
        mjd = None

        with open(inf_file, 'r') as file:
            lines = file.readlines()
            for line in lines:
                if epoch_key in line:
                    _, value = line.split('=', 1)
                    mjd = value.strip()
                    break

        t = Time(mjd, scale='utc', format='mjd')
        print(t.yday)

        # convert to required format
        input_format = "%Y:%j:%H:%M:%S.%f"
        dt = datetime.strptime(t.yday, input_format)

        # Format the datetime object using strftime
        output_format = "%y/%j_%H:%M:%S"
        formatted_time = dt.strftime(output_format)
        return t, formatted_time


def extract_pulses(sp_file):
    """
    Takes text file of pulses and returns list of time during 
    which pulses occured. Also return the binwidth of pulse

    spfile: single pulse file
    start_time: start time of data file
    """
    tt = []
    ww = []
    snr = []
    with open(sp_file) as fin:
        for line in fin:
            if line[0] in ["#", "\n"]:
                continue
            else: pass
            cols = line.split()
            
            tt.append(float(cols[2]))
            ww.append(int(cols[4]))
            snr.append(float(cols[1]))

    tt = np.array(tt)
    ww = np.array(ww)
    snr = np.array(snr)
            
    return tt, ww, snr


def write_offset(outfile, offsets, widths, snrs, times):
    """
    Write offsets to file
    """
    with open(outfile, 'w') as fout:
        hdr_fmt = "#{:^5s}  {:^12s}  {:^6s}  {:^8s}  {:^13s}\n"
        hdr = hdr_fmt.format("num", "offset", "w_bins", "snr", "time")  
        fout.write(hdr)
        fmt_str = "{:<5d}  {:>12.9f}  {:>6d}  {:>8.2f}  {:>13.6f}\n"
        for ii in range(len(offsets)):
            ostr = fmt_str.format(ii, offsets[ii], widths[ii], 
                                  snrs[ii], times[ii])
            fout.write(ostr)
    return 


def convert_times(times_after_start, start_time, record_time):
    """
    Takes pulse times in seconds and adds it to start_time 
    to return new list of times in correct format.
    ***We center the time such that the pulse occurs in the middle.
    Format of start_time is in YY/Day_HH:MM:SS
    """
    format_str = "%y/%j_%H:%M:%S"
    new_times = []

    for time in times_after_start:
        # convert to dt format
        dt = datetime.strptime(start_time, format_str)

        # add time and subtract half the record time to center it
        dt += timedelta(seconds=time)
        dt -= timedelta(seconds=record_time/2)
        new_times.append(dt.strftime(format_str))
    print(new_times)
    return new_times


def vrad_2_vdr(pulse_times, vrad_file, out_dir, 
               data_amount, srcdir=srcdir):
    """
    Converts raw telescope data (*.vrad) given list of times 
    of pulses to vdr format
    
    Executes script in format:

    rdef_io -V -f TIME -n 60 INPUT_FILE OUTPUT_FILE

    rdef_loc: location of rdef_io script
    """
    # Create directory to put vdr files in 
    vdr_dir = "%s/vdr" % (out_dir)
    if not os.path.exists(vdr_dir):
        os.makedirs(vdr_dir, exist_ok=True)

    for i, t in enumerate(pulse_times):
        # replace with vdr suffix to create output path + place 
        # in vdr folder
        basename = os.path.basename(vrad_file)
        basename = os.path.splitext(basename)[0]

        # zero pad pulse_str
        pulse_str = str(i).zfill(4)

        out_file = "%s/%s_p%s.vdr" % (vdr_dir, basename, pulse_str)

        # wait for process to finish before executing another
        rdef_loc = "%s/rdef_io" %srcdir
        cmd = "%s -V -f %s -n %d %s %s" % (rdef_loc, t, data_amount, 
                                           vrad_file, out_file)
        
        # if file exists already, we don't need to run
        if os.path.isfile(out_file):
            os.remove(out_file)
            # print("File: %s already exists" %out_file)
        else:
            print(cmd)
            subprocess.run(cmd, shell=True)


def vdr_2_cs(vdr_file, out_dir, freq_band, source_name, freq, bw, 
             telescope, srcdir=srcdir, band_num=1):
    """
    Converts *.vdr files to *.cs files
    Executes script in format:

    vsrdump_char_100 {vdr_file} -o {output_file} 
    -name {source_name} -comp -freq {freq} -bw {bw} -{telescope} -vsr
    """
    # Create directory to put cs files in 
    cs_dir = "%s/cs" % (out_dir)
    if not os.path.exists(cs_dir):
        os.makedirs(cs_dir, exist_ok=True)

    # replace with cs suffix to create output path + place in cs folder
    basename = os.path.basename(vdr_file)
    basename = os.path.splitext(basename)[0]

    # outfile naming convention will be xlcp_p0000-0001-01.cs
    pulse_num = re.search(r'p(\d{4})', vdr_file).group(0)
    polarization = ""
    if "rcp" in vdr_file.lower():
        polarization = "rcp"
    if "lcp" in vdr_file.lower():
        polarization = "lcp"

    out_file = "%s/cs/%s%s-%s-0%s.cs" % (out_dir, freq_band, 
                                      polarization, pulse_num, band_num)

    if freq_band == "s":
        vsrdump_loc = "%s/vsrdump_char_100" % srcdir
    else: 
        vsrdump_loc = "%s/vsrdump_char" % srcdir

    cmd = "%s %s -o %s -name %s -comp -freq %.2f -bw %.2f -%s -vsr" \
        % (vsrdump_loc, vdr_file, out_file, source_name, 
           freq, bw, telescope)

    # if file exists already, don't need to run
    if os.path.isfile(out_file):
        os.remove(out_file)
        # print("File: %s already exists" %out_file)
    else:
        subprocess.run(cmd, shell=True)

    return


def calculate_offsets(off_fname, mjd_start, times_after_start, 
                      bin_widths, snrs, cs_files, out_dir, 
                      freq_data, freq_sp, dm):
    """
    Calculates offset time for each burst.
    Writes offset to txt file for burst plotting
    
    delta_t = t_sp - t_cs
    """
    offsets = [0 for _ in range(len(times_after_start))]

    for i, time in enumerate(times_after_start):

        # calculate t_sp by adding pulse time to mjd start time
        dt = TimeDelta(time, format='sec')
        t_sp = mjd_start+dt

        # get cs file time
        hd, hsize, err = fb.read_header(cs_files[i], 4, fb.fmtdict)
        t_cs = Time(hd['tstart'], format='mjd')
        
        # subtract and return value in seconds
        off_set = t_sp-t_cs
        off_set = off_set.to_value('s')

        # Use regex to find the pattern. We are looking for p 
        # followed by 4 digits
        match = re.search(r'p(\d{4})', cs_files[i])
        index = -1
        if match:
            # Convert the digit part to an integer
            index = int(match.group(1))
        else:
            print("Pattern not found in filename.")
        
        # dispersive delay
        # t_dm = 4.15 * 10^3 s * (freq_lo^-2 – freq_hi^-2) * DM
        DMC = 4.148808e3

        t_dm = DMC * (freq_sp**-2. - freq_data**-2.) * dm

        # subtract from offset
        off_set -= t_dm

        offsets[index] = off_set

    off_file = "%s/%s" %(out_dir, off_fname)
    write_offset(off_file, offsets, bin_widths, snrs, times_after_start)
    return offsets


@click.command()
@click.option("--inf_file", type=str, 
              help="Full path of .inf file OR scan table", 
              required=True)
@click.option("--sp_file", type=str, 
              help="Full path of single pulse file", 
              required=True)
@click.option("--vrad_dir", type=str,
              help="Directory containing vrad files", 
              required=True)
@click.option("--vrad_base", type=str,
              help="Basename of vrad files {vrad_base}*.vrad", 
              required=True)
@click.option("--off_file", type=str, 
              help="Name of offset times file", required=True)
@click.option("--out_dir", type=str,
              help="Output directory of *.cs and intermediate *.vdr files",
              required=True)
@click.option("--freq_band", type=str,
              help="Frequency setup name of data", 
              required=True)
@click.option("--sp_freq_band", type=str,
              help="Frequency setup name of sp file" +\
                   "  (default: same as data)", 
              required=False)
@click.option("--dm", type=float,
              help="Dispersion Measure", 
              required = True)
@click.option("--source", type=str,
              help="Name of the pulse source", 
              required=True)
@click.option("--telescope", type=str,
              help="Name of telescope", 
              required=True)
@click.option("--data_amount", type=float,
              help="Amount of data to write (in seconds)",  
              required=True)
@click.option("--nproc", type=int,
              help="Number of processes to use in multiprocessing", 
              required=True)
def vrad_2_cs(inf_file, sp_file, off_file, vrad_dir, vrad_base, 
              out_dir, freq_band, dm, source, telescope, data_amount, 
              nproc, sp_freq_band=None, srcdir=srcdir):
    """
    Converts sections of *.vrad files into *.cs files around pulse times.

    Pulse times are located in a .singlepulse file. These times are equal to time after start

    Starting time is located in *.inf file. 
    The pulse times will have to be converted by adding the starting time.

    Input *.vrad files are in the form:
    
        {vrad_base}.{band}.vrad 
        Ex: 22-295-001_d63_PSR2_SLCP.vrad

    These *.vrad files will first be converted to *.vdr files and named:

        {out_dir/vdr/{base_name}_p{num}.vdr} 
        Ex: /data/vdr/22-295-001_d63_PSR2_SLCP_p1.vdr

    The *.vdr files will be converted to *.cs files and placed in

        {out_dir/cs/{burst_num}_{base_name}.cs}  
        Ex: /data/cs/22-295-001_d63_PSR2_SLCP_p1.cs

    """
    # dictionary containing freq for every freq band
    #freq_dct = {"s": 2250, "x1": 8224, "x2": 8256, "x3": 8288, "x4": 8320}
    #bw_dct = {"s": 100, "x": 32}

    # Get frequency setup for data and sp file
    fs_d = fs.get_freq_info(freq_band)
    if sp_freq_band is None:
        fs_sp = fs.get_freq_info(freq_band)
    else:
        fs_sp = fs.get_freq_info(sp_freq_band)

    # Extract starting time from .inf file
    mjd_start, start_time = extract_start_time(inf_file)
    times_after_start, bin_widths, snrs = extract_pulses(sp_file)

    # Add list of times to start time
    converted_times = convert_times(times_after_start, start_time, 
                                    data_amount)
    
    # Take each vrad file and convert to vdr
    print("\n\nvrad -> vdr...")
    vrad_infiles = []
    for pname in fs_d.pol_names:
        vgl = glob.glob("%s/%s*%s.vrad" %(vrad_dir, vrad_base, pname))
        vrad_infiles += vgl

    print(vrad_infiles)

    if len(vrad_infiles) == 0:
        print("No *vrad files found!")
        print("vrad_dir: %s" %vrad_dir)
        print("vrad_base: %s" %vrad_base)
        sys.exit(0)

    with mp.Pool(nproc) as pool:
        # Create a list of arguments for each call to vrad_2_vdr
        args_list = [(converted_times, vrad_file, out_dir, 
                      data_amount, srcdir) for vrad_file in vrad_infiles]

        # Use the pool to map the vrad_2_vdr function to the args_list
        pool.starmap(vrad_2_vdr, args_list)

    # Take vdr files and convert to cs
    print("\n\nvdr -> cs...")
    # We know all vdr files will be contained in vdr directory
    vdr_infiles = glob.glob("%s/vdr/*" % (out_dir))
    
    if len(vdr_infiles) == 0:
        print("No *vdr files found!")
        print("vdr_dir: %s/vdr" %out_dir)
        sys.exit(0)

    band_info = []
    bw_sb = fs_d.sub_bw 
    for vdr_file in vdr_infiles:
        band_num = re.search(r'-00(\d+)', vdr_file).group(1)
        bb_ii = int(band_num) - 1
        freq_sb = fs_d.sub_freqs[bb_ii]
        print(vdr_file, freq_sb, bw_sb, band_num)
        band_info.append((vdr_file, freq_sb, bw_sb, band_num))

    band = fs_d.band
    with mp.Pool(nproc) as pool:
        # Create a list of arguments for each call to vdr_2_cs
        args_list = [(vdr_file, out_dir, band, source, 
                    freq_sb, bw_sb, telescope, srcdir, band_num) 
                    for vdr_file, freq_sb, bw_sb, band_num in band_info]

        # Use the pool to map the vdr_2_cs function to the args_list
        pool.starmap(vdr_2_cs, args_list)

    # find cs files to calculate offset times
    cs_infiles = sorted(glob.glob("%s/cs/*" % (out_dir)))
    print(times_after_start)
    print(cs_infiles)
    
    if len(cs_infiles) == 0:
        print("No *cs files found!")
        print("cs_dir: %s/cs" %out_dir)
        sys.exit(0)

    # take subset of cs_infiles for each pulse
    # not the most efficient way of doing it
    subset = []
    for i in range(len(times_after_start)):
        for cfile in cs_infiles:
            if int(re.search(r'p(\d{4})', cfile).group(1)) == i:
                subset.append(cfile)
                break
    cs_infiles = subset
    print(cs_infiles)

    # Get top of band freq for data and for sp file
    fhi_d  = fs_d.freq_hi()
    fhi_sp = fs_sp.freq_hi()
    
    calculate_offsets(off_file, mjd_start, times_after_start, 
                      bin_widths, snrs, cs_infiles, out_dir, 
                      fhi_d, fhi_sp, dm)

    return 

debug = 0

if __name__ == "__main__":
    if not debug: 
        vrad_2_cs()
