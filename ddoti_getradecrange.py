#!/usr/bin/python
"""
 ddoti_getradecrange file_list [weight_min]
"""

from pyfits import getdata
from numpy import where,loadtxt
from fit_wcs import xy2ad
import sys,os

def usage():
    print __doc__
    sys.exit()

def ddoti_getradecrange(file_list,wt_min=0.):
    """
    work out the ra and dec range covered by each file
    """

    filename, file_extension = os.path.splitext(file_list)
    if (file_extension=='.fits'): files=[file_list]
    else: files = loadtxt(file_list,dtype='string').tolist()

    for file in files:

        wfile=file.replace('.fits','.wt.fits')
        if (os.path.exists(wfile)==False): usage()

        wt,hdr = getdata(wfile,header=True)
        wtm = wt.max()
        wtx,wty = wt.mean(axis=0), wt.mean(axis=1)

        wx,wy = where(wtx>wt_min*wtm)[0],where(wty>wt_min*wtm)[0]

        x1,x2 = wx[0]+1,wx[-1]+1
        y1,y2 = wy[0]+1,wy[-1]+1

        r1,d1 = xy2ad(x1,y1,hdr)
        r2,d2 = xy2ad(x2,y2,hdr)

        if (r2<r1): r1,r2=r2,r1
        if (d2<d1): d1,d2=d2,d1

        print """sethead I1=%d I2=%d J1=%d J2=%d RA0=%.8f RA1=%.8f DEC0=%.8f DEC1=%.8f %s""" % (x1,x2,y1,y2,r1,r2,d1,d2,file)


if __name__ == '__main__':
    """
    """
    if (len(sys.argv)<2): usage()

    file_list = sys.argv[1]
    if (os.path.exists(file_list)==False): usage()

    wt_min=0.
    if (len(sys.argv)>2): wt_min=float(sys.argv[2])

    ddoti_getradecrange(file_list,wt_min=wt_min)
