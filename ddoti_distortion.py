#!/usr/bin/python
"""
 ddoti_distortion.py file_list
"""
import sys,os

import pyfits
from scipy.optimize import fmin
from fit_wcs import ij2ad,ij2xy,xy2ij
from linfit import linfit
from numpy import loadtxt,cos,array,sin,sqrt,eye,dot,vstack,zeros

import matplotlib as mpl
mpl.use('Agg')

from matplotlib.pyplot import plot,xlabel,ylabel,savefig,legend,grid,title

def usage():
    print __doc__
    sys.exit()

def ddoti_distortion(file_list,nmax=1000):
    """
    """
    fac = 57.29577951

    filename, file_extension = os.path.splitext(file_list)
    if (file_extension=='.fits'): files=[file_list]
    else: files = loadtxt(file_list,dtype='string',ndmin=1).tolist()

    nf0 = len(files)*nmax
    r0 = zeros(nf0,dtype='float64')
    d0 = zeros(nf0,dtype='float64')
    r = zeros(nf0,dtype='float64')
    d = zeros(nf0,dtype='float64')
    ra0 = zeros(nf0,dtype='float64')
    dec0 = zeros(nf0,dtype='float64')
    c = zeros(nf0,dtype='float64')
    i = zeros(nf0,dtype='float64')
    j = zeros(nf0,dtype='float64')

    k0=0
    for file in files:
        base=file.replace('.fits','')
        srcfile=base+'_dir/'+base+'_radec.txt.match.txt'
        dat = loadtxt(srcfile,usecols=(2,4,5,8,9)).T
        ii = dat[0].argsort()[:nmax]
        k1 = k0 + len(ii)
        r0[k0:k1],d0[k0:k1],x,y = dat[1:,ii]
        hdr = pyfits.getheader(file)
        i[k0:k1],j[k0:k1] = xy2ij(x,y,hdr,distort=False)
        r[k0:k1],d[k0:k1] = ij2ad(i[k0:k1],j[k0:k1],hdr)
        ra0[k0:k1],dec0[k0:k1] = hdr['CRVAL1'],hdr['CRVAL2']
        c[k0:k1]=cos(hdr['CRVAL2']/fac)
        k0+=len(ii)

    if (k0<nf0):
        r0 = r0[:k0]; d0 = d0[:k0]; r = r[:k0]; d = d[:k0]; ra0 = ra0[:k0]
        dec0 = dec0[:k0]; c = c[:k0]; i = i[:k0]; j = j[:k0]

    try: a0,b0 = hdr['PV1_7']*3600.,hdr['PV1_17']*3600.
    except: a0,b0 = 1.4,0.1
    res = array([0.,0.,a0,b0])

    sd0,cd0 = sin(d0/fac),cos(d0/fac)
    def myfun(par,i,j):
        i00,j00 = par[0]/360., par[1]/360.
        #
        aa0, dd0 = (i00/c+ra0)/fac, (j00+dec0)/fac
        sdd0,cdd0,cfac = sin(dd0),cos(dd0),cos(r0/fac-aa0)
        cos_c = sdd0*sd0 + cdd0*cd0*cfac
        i0 = fac*cd0*sin(r0/fac-aa0) / cos_c
        j0 = fac*(cdd0*sd0-sdd0*cd0*cfac) / cos_c
        #
        rad = sqrt( (i-i00)**2+(j-j00)**2 )
        rad0 = sqrt( i0**2+j0**2 )
        dis = rad0-rad
        dism = rad**3*(par[2] + par[3]*rad**2)/3600.
        return ( (dis-dism)**2 ).mean()


    matr = eye(2)
    for niter in xrange(5):
        par0 = res
        res = fmin(myfun,par0,args=(i,j),disp=False); res = fmin(myfun,res,args=(i,j),disp=False)
        #
        i00,j00 = res[0]/360.,res[1]/360.
        aa0, dd0 = (i00/c+ra0)/fac, (j00+dec0)/fac
        sdd0,cdd0,cfac = sin(dd0),cos(dd0),cos(r0/fac-aa0)
        cos_c = sdd0*sd0 + cdd0*cd0*cfac
        i0 = fac*cd0*sin(r0/fac-aa0) / cos_c
        j0 = fac*(cdd0*sd0-sdd0*cd0*cfac) / cos_c
        #
        rad2 = (i-i00)**2+(j-j00)**2
        i1 = i + (i-i00)*rad2*(res[2]+res[3]*rad2)/3600.
        j1 = j + (j-j00)*rad2*(res[2]+res[3]*rad2)/3600.
        #
        rii,rji = linfit(i1-i00,i1-i00-i0), linfit(i1-i00,j1-j00-j0)
        rjj,rij = linfit(j1-j00,j1-j00-j0), linfit(j1-j00,i1-i00-i0)
        rsd1,rsd2 = abs(rji[1])*3600.,abs(rij[1])*3600.
        print """ Iter %d %.8f %.8f Resids: %.8f %.8f (arcsec) %s""" % (niter,res[2]/3600.,res[3]/3600.,rsd1,rsd2,files[0])
        #
        matr0 = array([[1-rii[1],-rij[1]],[-rji[1],1-rjj[1]]])
        matr = dot(matr0,matr)
        i,j = dot(matr0,vstack((i,j)))
        if (rsd1<1.e-3 and rsd2<1.e-3): break


    # diagnostic plot
    rad2 = (i-i00)**2+(j-j00)**2
    i1 = i + (i-i00)*rad2*(res[2]+res[3]*rad2)/3600.
    j1 = j + (j-j00)*rad2*(res[2]+res[3]*rad2)/3600.
    rad20 = i0**2+j0**2
    dis=sqrt(rad20)-sqrt(rad2)
    dis1=sqrt(rad20)-sqrt((i1-i00)**2+(j1-j00)**2)
    plot (sqrt(rad2),dis*3600.,'o',label='Uncorrected',alpha=0.3); plot (sqrt(rad2),dis1*3600.,'o',label='Corrected',alpha=0.3)
    grid()
    title("""Distortion Params: %.8f %.8f\n""" % (res[2]/3600.,res[3]/3600.),fontsize=16)
    xlabel("Radius from Image Center [Degrees]",fontsize=16)
    ylabel("Source Offset [arcsec]",fontsize=16)
    legend(loc=2)
    savefig(files[0].replace('.fits','_dist.jpg'))

    # a*rad^3 + b*rad^5
    for fitsfile in files:
        a,b = res[2]/3600.,res[3]/3600.
        os.system("""sethead PV1_1=1.0 PV2_1=1.0 %s""" % fitsfile)
        os.system("""sethead PV1_7=%.8f PV1_9=%.8f PV2_7=%.8f PV2_9=%.8f %s""" % (a,a,a,a,fitsfile))
        os.system("""sethead PV1_17=%.8f PV1_21=%.8f PV2_17=%.8f PV2_21=%.8f %s""" % (b,b,b,b,fitsfile))
        os.system("""sethead PV1_19=%.8f PV2_19=%.8f %s""" % (2*b,2*b,fitsfile))
        hdr = pyfits.getheader(fitsfile)
        cd = array([[hdr['CD1_1'],hdr['CD1_2']],[hdr['CD2_1'],hdr['CD2_2']]])
        cd1 = dot(matr,cd)
        os.system("""sethead CD1_1=%.8f CD1_2=%.8f CD2_1=%.8f CD2_2=%.8f %s""" % (cd1[0,0],cd1[0,1],cd1[1,0],cd1[1,1],fitsfile))
        hdr = pyfits.getheader(fitsfile)
        x,y = ij2xy(i00,j00,hdr,distort=True)
        r,d = ij2ad(i00,j00,hdr)
        os.system("""sethead CRPIX1=%.2f CRPIX2=%.2f CRVAL1=%.8f CRVAL2=%.8f %s""" % (x,y,r,d,fitsfile))


if __name__ == "__main__":
    """
    """
    if (len(sys.argv)<2): usage()

    file_list=sys.argv[1]
    if (os.path.exists(file_list)==False): usage()

    ddoti_distortion(file_list)
