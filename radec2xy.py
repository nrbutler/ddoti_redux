#!/usr/bin/python
"""
radec2xy.py fits_file radec_file
"""

import sys,os
from fit_wcs import ad2xy
import pyfits
from numpy import array

def usage():
    print __doc__
    sys.exit()

def radec2xy(fits_file,radecfile):
    """
    """
    hdr=pyfits.getheader(fits_file)
    distort = 'PV1_1' in hdr

    f = open(radecfile,'r')
    ra=[];dec=[]
    str=[]
    for line in f.readlines():
        str.append(line.strip())
        ll=line.split()
        if (len(ll)<=1): continue
        if (ll[0][0]=="#"):
            ra.append(-999)
            dec.append(-999)
        else:
            ra.append(eval(ll[0]))
            dec.append(eval(ll[1]))

    f.close()
    ra = array(ra); dec = array(dec)
    x,y = 1.*ra,1.*dec
    j = (ra!=-999)*(dec!=-999)
    x[j],y[j] = ad2xy(ra[j],dec[j],hdr,distort=distort)

    for i in xrange(len(x)):
        if (ra[i]!=-999 and dec[i]!=-999):
            print """%f %f""" % (x[i],y[i])
        else: print str[i]

def main():

    if (len(sys.argv)<3): usage()

    fits_file=sys.argv[1]
    radecfile=sys.argv[2]
    radec2xy(fits_file,radecfile)

if __name__ == "__main__":
    main()
