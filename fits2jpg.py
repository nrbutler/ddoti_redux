#!/usr/bin/python
"""
 fits2jpg.py fitsfile [size] [noinvert]
"""

import sys,os

from scipy.misc import imresize,imsave

import pyfits
from numpy.random import rand
from numpy import sqrt

def usage():
    print __doc__
    sys.exit()

def fits2jpg(file0,nsample=50000,out_size=1024,wt_min=1.e-3,invert=True):
    """
    """
    file1=file0.replace('.fits','.jpg')

    outfile=file0.replace('.fits','.jpg')
    wfile=file0.replace('.fits','.wt.fits')
    sys.stderr.write("""Creating %s from %s ...\n""" % (outfile,file0))

    x=pyfits.getdata(file0)
    a,b=x.shape
    n=a*b

    if (os.path.exists(wfile)):
        wt = pyfits.getdata(wfile)
        h = wt>wt_min
        x[~h]=0
        x[h]*=sqrt(wt[h])

    if (nsample>n): nsample=n
    xind = (rand(nsample)*(a-1)+1).astype('int32')
    yind = (rand(nsample)*(b-1)+1).astype('int32')

    # summary jpg images
    x1 = x[xind,yind]
    s= x1.argsort()
    s1,s2 = s[int(0.25*nsample)],s[int(0.95*nsample)]
    x1,x2 = x1[s1],x1[s2]

    if (invert): y = 255.*(x1-x.clip(x1,x2))/(x2-x1)
    else: y = 255.*(x.clip(x1,x2)-x1)/(x2-x1)
    if (b!=out_size):
        nx,ny = int(round(1.*out_size*a/b)),out_size
        sys.stderr.write("""Resizing Image %s to %dx%d\n""" % (file1,ny,nx))
        y = imresize(y,size=(nx,ny),mode='F')

    imsave(file1,y[::-1])


if __name__ == "__main__":
    """
    """
    if (len(sys.argv)<2): usage()

    file0=sys.argv[1]
    if (os.path.exists(file0)==False): usage()

    sz=1024
    if (len(sys.argv)>2): sz=int(sys.argv[2])

    invert=True
    if (len(sys.argv)>3):
        if(sys.argv[3]=="noinvert"): invert=False

    fits2jpg(file0,out_size=sz,invert=invert)
