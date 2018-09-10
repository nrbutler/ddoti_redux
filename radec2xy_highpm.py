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
    if (len(data)>=6): dra,ddec = data[4:6]
    else: dra,ddec = 0*ra,0*ra

    x,y = ad2xy(ra,dec,hdr,distort=distort)
    nx=hdr['NAXIS1']; ny=hdr['NAXIS2']
    h = (x>=1)*(x<=nx)*(y>=1)*(y<=ny)
    x = x[h]; y = y[h]; mag = mag[h]; dmag=dmag[h]; ra=ra[h]; dec=dec[h]; dra=dra[h]; ddec=ddec[h]
    if (len(data)>=6):
        x1,y1 = ad2xy(ra+dra,dec,hdr,distort=distort)
        x2,y2 = ad2xy(ra,dec+ddec,hdr,distort=distort)
        for i in xrange(len(x)): print """%f %f %d %f %f %f %f""" % (x[i],y[i],i+1,mag[i],dmag[i],x1[i]+x2[i]-2*x[i],y1[i]+y2[i]-2*y[i])
    else:
        for i in xrange(len(x)): print """%f %f %d %f %f""" % (x[i],y[i],i+1,mag[i],dmag[i])


    # consider a high proper motion star with path x,y = x0+a*dx,y0+a*dy
    #   find all stars within some decently sized radius from x0,y0, at positions x1,y1
    # the shortest distance from each new point x1,y1 to propoer motion star's path is: dis = [dy*(x1-x0)-dx*(y1-y0)]/sqrt( dx^2+dy^2 )

    #, for each star x1,y1 which we want to compare to a set of high proper motion stars at x0,y0, calculate dis1 = [dy*(x1-x0)-dx*(y1-y0)]/sqrt( dx^2+dy^2 )
    #    if that's below some threshold AND dis1_0 = sqrt( (x1-x0)^2+(y1-y0)^2 ) is also below some larger threshold, then nuke it
    #   perhaps dis1_thresh = 5 pixels and dis1_0_thresh = 10 pixels


def main():

    if (len(sys.argv)<3): usage()

    fits_file=sys.argv[1]
    radecfile=sys.argv[2]
    radec2xy(fits_file,radecfile)

if __name__ == "__main__":
    main()
