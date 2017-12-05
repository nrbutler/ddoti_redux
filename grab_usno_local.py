#!/usr/bin/python
"""
grab_usno_local.py fits_file ra dra dec ddec
"""

import sys,os
import pyfits
try: from weave import inline
except: from scipy.weave import inline
from numpy import array,arange

def usage():
    print __doc__
    sys.exit()


def grab_usno_local(fitsfile,ra1,ra2,dec1,dec2):
    """
    """
    dat = pyfits.getdata(fitsfile).astype('float64')
    r,d,m = dat[0],dat[1],dat[2]

    rdat = array([ra1,ra2],dtype='float64')
    idat = array([0,0],dtype='int64')
    n=len(r)

    search_code = """
    unsigned long i0,i1;
    for (i0=0;i0<n;i0++) {
        if (r[i0]>=rdat[0]) break;
    }
    idat[0] = i0;
    for (i1=i0;i1<n;i1++) {
        if (r[i1]>=rdat[1]) break;
    }
    idat[1] = i1-1;
    """
    #
    # find matches at a resolution of 2/fwhm/sqrt(2)
    inline(search_code,['n','r','rdat','idat'])
    i0,i1 = idat
    r,d,m = r[i0:i1+1],d[i0:i1+1],m[i0:i1+1]

    h = (d>=dec1)*(d<=dec2)
    r,d,m = r[h],d[h],m[h]
    dm = 0*m+999
    h = m<=21
    dm[h] = 0.05 * 10**(0.11*(m[h]-10).clip(0))
    for i in xrange(len(r)):
        print """%.8f %.8f %.4f %.4f""" % (r[i],d[i],m[i],dm[i])


def main():
    """
    """
    if (len(sys.argv)<5): usage()

    fits_file=sys.argv[1]
    if (os.path.exists(fits_file)==0): usage()

    file0 = os.path.basename(fits_file)
    fs = file0.split('_')
    ra0,dec0 = float(fs[2]), float(fs[3].strip('.fits'))

    ra=float(sys.argv[2])
    dra=float(sys.argv[3])
    dec=float(sys.argv[4])
    ddec=float(sys.argv[5])

    ra1,ra2 = ra-dra/2.,ra+dra/2.
    if (ra1<0 and ra0>180): ra1,ra2 = ra1+360,ra2+360

    dec1,dec2 = dec-ddec/2.,dec+ddec/2.
    if (dec1<-90 and dec0>0): dec1,dec2 = dec1+180,dec2+180

    grab_usno_local(fits_file,ra1,ra2,dec1,dec2)


if __name__ == "__main__":
    main()
