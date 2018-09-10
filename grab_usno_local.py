#!/usr/bin/python
"""
grab_usno_local.py fits_file outfile ra dra dec ddec [year]
"""

import sys,os
from pyfits import getdata,writeto
try: from weave import inline
except: from scipy.weave import inline
from numpy import array,vstack,cos,zeros

def usage():
    print __doc__
    sys.exit()


def grab_usno_local(fitsfile,outfile,ra1,ra2,dec1,dec2,year=2018.):
    """
    """
    r,d,m = getdata(fitsfile).astype('float32')

    prop_motion=False
    fitsfile1 = fitsfile.replace('.fits','.pm.fits')
    if (os.path.exists(fitsfile1)):
        prop_motion=True
        dr,dd,pm = 0.*r,0.*d,zeros(len(r),dtype='int16')
        sys.stderr.write("""Making proper motion corrections for %s.\n""" % fitsfile)
        ii,dx,dy = getdata(fitsfile1).astype('int32')
        pm[ii]=1
        cd=cos(0.5*(dec1+dec2)/57.296)
        dr[ii] = dx*(year-2000.)/3.6e6
        dd[ii] = dy*(year-2000.)/3.6e6
        r[ii] += dr[ii]/cd
        d[ii] += dd[ii]

    rdat = array([ra1,ra2],dtype='float32')
    idat = array([0,0],dtype='int32')
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
    if (prop_motion): dr,dd,pm = dr[i0:i1+1],dd[i0:i1+1],pm[i0:i1+1]

    h = (d>=dec1)*(d<=dec2)*(r>=ra1)*(r<=ra2)
    r,d,m = r[h],d[h],m[h]
    if (prop_motion): dr,dd,pm = dr[h],dd[h],pm[h]
    dm = 0*m+999
    h = m<=21
    dm[h] = 0.05 * 10**(0.11*(m[h]-10).clip(0))

    writeto(outfile, vstack((r,d,m,dm)) )
    if (prop_motion):
        outfile1 = outfile.replace('.fits','.pm.fits')
        h = pm==1
        if (h.sum()>0): writeto(outfile1, vstack((r[h],d[h],m[h],dm[h],dr[h],dd[h])) )


def main():
    """
    """
    if (len(sys.argv)<7): usage()

    fits_file=sys.argv[1]
    if (os.path.exists(fits_file)==0): usage()

    outfile=sys.argv[2]

    file0 = os.path.basename(fits_file)
    fs = file0.split('_')
    ra0,dec0 = float(fs[2]), float(fs[3].strip('.fits'))

    ra=float(sys.argv[3])
    dra=float(sys.argv[4])
    dec=float(sys.argv[5])
    ddec=float(sys.argv[6])

    year=2018.
    if (len(sys.argv)>7): year=float(sys.argv[7])

    ra1,ra2 = ra-dra/2.,ra+dra/2.
    if (ra1<0 and ra0>180): ra1,ra2 = ra1+360,ra2+360

    dec1,dec2 = dec-ddec/2.,dec+ddec/2.
    if (dec1<-90 and dec0>0): dec1,dec2 = dec1+180,dec2+180

    grab_usno_local(fits_file,outfile,ra1,ra2,dec1,dec2,year=year)


if __name__ == "__main__":
    main()
