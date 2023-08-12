# Here is where we specify frequency setups

class FreqSetup():
    """
    Class object to hold freq info
    """
    def __init__(self):
        self.name = ""
        self.band = ""
        self.pol_names = []
        self.sub_freqs = []
        self.sub_bw = 0

    def freq_hi(self):
        if self.sub_bw > 0:
            fhi = max(self.sub_freqs) + 0.5 * self.sub_bw
        else:
            fhi = max(self.sub_freqs) - 0.5 * self.sub_bw
        return fhi
    
    def freq_lo(self):
        if self.sub_bw > 0:
            flo = min(self.sub_freqs) - 0.5 * self.sub_bw
        else:
            flo = min(self.sub_freqs) + 0.5 * self.sub_bw
        return flo

    def freq_mid(self):
        flo = self.freq_lo()
        fhi = self.freq_hi()
        fmid = 0.5 * (flo + fhi)
        return fmid

    def __repr__(self):
        return "FreqSetup(%s)" %self.name

    def __str__(self):
        return "FreqSetup(%s)" %self.name


def get_freq_info(fkey):
    """
    `fkey` is the name of the frequency setup 
    you want to use from those specified below.

    freq_dct gives a dictionary specifying the 
    center frequencies of the subbands

    sub_bw sets the bandwidth of the sub-band
    """
    # to avoid case issues make it lower 
    fkey = fkey.lower()

    valid_keys = ["sband-1", "xband-1"]

    fsetup = FreqSetup()
    
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

    else:
        print("Invalid key: %s" %fkey)
        print("Must be in %s" %(",".join(valid_keys)))
        return 0

    return fsetup
    
