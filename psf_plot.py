#!/usr/bin/python
"""
 psf_plot.py psf_string [cam]
"""

import sys,os

from scipy.interpolate import RectBivariateSpline
import pyfits
from glob import glob
from numpy import linspace,zeros,sqrt

import matplotlib as mpl
mpl.use('Agg')

from matplotlib.pyplot import imshow,xlabel,ylabel,savefig,title,grid,annotate

def psf_plot(psf_string='f\*stack\*dir/psf.fits',cam='C0'):
    """
    """
    psfs=glob(psf_string)
    psfs.sort()

    n=len(psfs)
    sz=zeros(n,dtype='int16')
    fwhm=zeros(n,dtype='float32')
    osamp=zeros(n,dtype='int16')

    psf=[]
    for i in xrange(n):
        d,h=pyfits.getdata(psfs[i],header=True)
        fwhm[i]=h['FWHM']
        sz[i]=h['NAXIS1']
        osamp[i]=h['OVERSAM']
        psf.append(d)

    sz0 = (sz-1)/(2*osamp)
    sz0m = sz0.max()
    mfwhm = 2*int(round(fwhm.max()))
    xx = linspace(-sz0m+mfwhm,sz0m-mfwhm,2*(sz0m-mfwhm)+1)
    nx = len(xx)
    ns = int(sqrt(n))
    sx = nx*ns
    psfa = zeros((sx,sx),dtype='float32')
    for i in xrange(n):
        i0 = int(psfs[i][1])
        j0 = int(psfs[i][2])
        xx0 = linspace(-sz0[i],sz0[i],sz[i])
        psf_s = RectBivariateSpline(xx0,xx0,psf[i])(xx,xx)
        psfa[i0*nx:(i0+1)*nx,j0*nx:(j0+1)*nx] = psf_s/psf_s.max()


    imshow(psfa,cmap='gray_r',origin='lower',interpolation='None')
    xlabel("Image X",fontsize=16); ylabel("Image y",fontsize=16)
    title(cam+" PSF Summary"); grid()

    for i in xrange(n):
        i0 = int(psfs[i][1])
        j0 = int(psfs[i][2])
        annotate("""%.1f""" % fwhm[i],[j0*nx,i0*nx],fontsize=16,alpha=0.3)

    savefig(cam+'_'+"psf_plot.jpg")


if __name__ == "__main__":
    """
    """
    psf_string=sys.argv[1]
    cam='C0'
    if (len(sys.argv)>2): cam=sys.argv[2]

    psf_plot(psf_string,cam=cam)
