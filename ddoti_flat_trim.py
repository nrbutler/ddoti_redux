#!/usr/bin/python
"""
 ddoti_flat_trim.py fits_file
"""
import sys,os

import pyfits
from scipy.signal import medfilt2d
from numpy import median,abs,logical_or

def usage():
    print __doc__
    sys.exit()

def ddoti_flat_trim(fits_file,nsigma=5,flat_min=0.3):
    """
    """
    dat = pyfits.getdata(fits_file)
    hdr = pyfits.getheader(fits_file)
    dat1 = dat - medfilt2d(dat,5)

    d0 = median(dat1)
    dat1 -= d0
    dat1 = abs(dat1)
    dd0 = 1.48*median(dat1)

    h = dat1>nsigma*dd0
    sys.stderr.write("""Flagging %d bad pixels in %s\n""" % (h.sum(),fits_file))
    dat[h]=0
    dat[dat<flat_min]=0

    os.remove(fits_file)
    hdr['FLATMIN'] = flat_min
    pyfits.writeto(fits_file,dat,hdr)


if __name__ == "__main__":
    """
    """
    if (len(sys.argv)<2): usage()

    fits_file=sys.argv[1]
    if (os.path.exists(fits_file)==False): usage()

    nsigma=5
    if (len(sys.argv)>2): nsigma=int(sys.argv[2])

    ddoti_flat_trim(fits_file,nsigma=nsigma)
