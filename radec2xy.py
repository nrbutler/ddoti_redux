#!/usr/bin/python
"""
radec2xy.py fits_file radec_file
"""

import sys,os
from fit_wcs import ad2xy
from pyfits import getheader,getdata
from numpy import loadtxt

def usage():
    print __doc__
    sys.exit()

def radec2xy(fits_file,radecfile):
    """
    """
    hdr=getheader(fits_file)
    distort = 'PV1_1' in hdr

    filename, file_extension = os.path.splitext(radecfile)
    if (file_extension=='.fits'): data = getdata(radecfile)
    else: data = loadtxt(radecfile,ndmin=2).T

    ra,dec = data[:2]
    if (len(data)>=4): mag,dmag = data[2:4]
    else: mag,dmag = 0*ra+19,0*dec+0.1

    x,y = ad2xy(ra,dec,hdr,distort=distort)
    nx=hdr['NAXIS1']; ny=hdr['NAXIS2']
    h = (x>=1)*(x<=nx)*(y>=1)*(y<=ny)
    x = x[h]; y = y[h]; mag = mag[h]; dmag=dmag[h]; ra=ra[h]; dec=dec[h]
    for i in xrange(len(x)): print """%f %f %d %f %f""" % (x[i],y[i],i+1,mag[i],dmag[i])


def main():

    if (len(sys.argv)<3): usage()

    fits_file=sys.argv[1]
    radecfile=sys.argv[2]
    radec2xy(fits_file,radecfile)

if __name__ == "__main__":
    main()
