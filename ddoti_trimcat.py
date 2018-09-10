#!/usr/bin/python
"""
ddoti_trimcat.py phot_file cat_file [nmax] [snr] [phot_rad] [cat_rad] [bag_bright]
"""

import sys,os
from numpy import loadtxt,ones,cos,median
from pyfits import getdata

def usage():
    print __doc__
    sys.exit()


def ddoti_trimcat(phot_file,cat_file,phot_rad=30.0,cat_rad=90.0,snr=10.,nmax=100,mag_bright=6.):
    """
      only retain catalog sources with >snr, and
        not within phot_rad of another
        not within cat_rad of another bright catalogued source
    """

    phot_dat = loadtxt(phot_file,ndmin=2).T
    if (len(phot_dat[0])==0): sys.exit()

    cat_dat = getdata(cat_file)
    r0,d0,m0 = cat_dat[:3,cat_dat[2]<mag_bright+6]

    phot_dat = phot_dat[:,phot_dat[3]<1./snr]
    if (len(phot_dat[0])>0):

        min_rad2 = (0.1/3600)**2
        phot_rad2 = (phot_rad/3600.)**2
        if (phot_rad2<min_rad2): phot_rad2 = 2*min_rad2
        cat_rad2 = (cat_rad/3600.)**2
        if (cat_rad2<min_rad2): cat_rad2 = 2*min_rad2

        cat_rad2 *= 10**(-0.3*(m0-mag_bright))

        r1,d1 = phot_dat[:2]
        cd = cos(median(d1)/57.3)
        n = len(r1)
        good = ones(n,dtype='bool')
        if (n>1):
            for i in xrange(n):
                dis2 = ((r1-r1[i])*cd)**2 + (d1-d1[i])**2
                dis2m = dis2[dis2>min_rad2].min()
                if (dis2m<phot_rad2): good[i] = False

            phot_dat = phot_dat[:,good]
            ii = phot_dat[2].argsort()
            phot_dat = phot_dat[:,ii]

        if (len(r0)>0):
            r1,d1 = phot_dat[:2]
            n = len(r1)
            good = ones(n,dtype='bool')
            if (n>1):
                for i in xrange(n):
                    dis2 = ((r0-r1[i])*cd)**2 + (d0-d1[i])**2
                    h = dis2>min_rad2
                    i0 = (dis2[h]-cat_rad2[h]).argmin()
                    dis2m = dis2[h][i0]
                    if (dis2m<cat_rad2[h][i0]): good[i] = False

                phot_dat = phot_dat[:,good]

        of=open(phot_file,'r')
        ln1 = of.readline()
        of.close()

        of=open(phot_file,'w')
        if (ln1[0]=='#'): of.write(ln1)
        n0 = len(phot_dat[0])
        for i in xrange(n0):
            if (i==nmax):
                print """Warning: only displaying at most the brightest %d of %d sources!""" % (nmax,n0)
                break
            ra,dec,mag,dmag,mag_big,dmag_big,fwhm,x,y,xa,ya,x2a,y2a,expos,idx = phot_dat[:,i]
            of.write("""%f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n""" % (ra,dec,mag,dmag,mag_big,dmag_big,fwhm,x,y,xa,ya,x2a,y2a,expos,i+1))

        of.close()


if __name__ == "__main__":
    """
    """

    if (len(sys.argv)<3): usage()

    phot_file=sys.argv[1]
    if (os.path.exists(phot_file)==0): usage()

    cat_file=sys.argv[2]
    if (os.path.exists(cat_file)==0): usage()

    nmax=100
    if (len(sys.argv)>3): nmax = int(sys.argv[3])
    snr=10.
    if (len(sys.argv)>4): snr = float(sys.argv[4])
    phot_rad=30.
    if (len(sys.argv)>5): phot_rad = float(sys.argv[5])
    cat_rad=90.
    if (len(sys.argv)>6): cat_rad = float(sys.argv[6])
    mag_bright=6.
    if (len(sys.argv)>7): mag_bright = float(sys.argv[7])

    ddoti_trimcat(phot_file,cat_file,phot_rad=phot_rad,cat_rad=cat_rad,snr=snr,nmax=nmax,mag_bright=mag_bright)
