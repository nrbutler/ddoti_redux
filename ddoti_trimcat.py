#!/usr/bin/python
"""
ddoti_trimcat.py cat_file file [nmax] [snr] [cat_rad]
"""

import sys,os
from numpy import loadtxt,ones,cos,median

def usage():
    print __doc__
    sys.exit()


def ddoti_trimcat(cat_file,cat_rad=30.0,snr=10.,nmax=100):
    """
      only retain catalog sources with >snr, and
        not within cat_rad of another
    """

    phot_dat = loadtxt(cat_file,ndmin=2).T
    if (len(phot_dat[0])==0): sys.exit()

    phot_dat = phot_dat[:,phot_dat[3]<1./snr]
    if (len(phot_dat[0])>0):
        r1,d1 = phot_dat[:2]

        n = len(r1)
        good = ones(n,dtype='bool')
        cd = cos(median(d1)/57.3)

        min_rad2 = (0.1/3600)**2
        cat_rad2 = (cat_rad/3600.)**2
        if (cat_rad2<min_rad2): cat_rad2 = 2*min_rad2

        if (n>1):
            for i in xrange(n):
                dis2 = ((r1-r1[i])*cd)**2 + (d1-d1[i])**2
                dis2m = dis2[dis2>min_rad2].min()
                if (dis2m<cat_rad2): good[i] = False

            phot_dat = phot_dat[:,good]
            ii = phot_dat[2].argsort()
            phot_dat = phot_dat[:,ii]

        of=open(cat_file,'r')
        ln1 = of.readline()
        of.close()

        of=open(cat_file,'w')
        if (ln1[0]=='#'): of.write(ln1)
        n0 = len(phot_dat[0])
        for i in xrange(n0):
            if (i==nmax):
                print """Warning: only displaying the brightest %d of %d sources!""" % (nmax,n0)
                break
            ra,dec,mag,dmag,mag_big,dmag_big,fwhm,x,y,xa,ya,x2a,y2a,expos,idx = phot_dat[:,i]
            of.write("""%f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n""" % (ra,dec,mag,dmag,mag_big,dmag_big,fwhm,x,y,xa,ya,x2a,y2a,expos,i+1))

        of.close()


if __name__ == "__main__":
    """
    """

    if (len(sys.argv)<2): usage()

    cat_file=sys.argv[1]
    if (os.path.exists(cat_file)==0): usage()

    nmax=100
    if (len(sys.argv)>2): nmax = int(sys.argv[2])
    snr=10.
    if (len(sys.argv)>3): snr = float(sys.argv[3])
    cat_rad=30.
    if (len(sys.argv)>4): cat_rad = float(sys.argv[4])

    ddoti_trimcat(cat_file,cat_rad=cat_rad,snr=snr,nmax=nmax)
