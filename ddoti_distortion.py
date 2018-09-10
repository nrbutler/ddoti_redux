#!/usr/bin/python
"""
 ddoti_distortion.py file_list [p1]
"""
import sys,os

from pyfits import getheader
from scipy.optimize import fmin
from fit_wcs import xy2ad,ij2xy,xy2ij
from linfit import linfit
from numpy import loadtxt,sin,cos,array,sqrt,dot
from nsigma_clip import nsigma_clip

import matplotlib as mpl
mpl.use('Agg')

from matplotlib.pyplot import plot,xlabel,ylabel,savefig,legend,grid,ylim

def usage():
    print __doc__
    sys.exit()

def ddoti_distortion(fits_file,nmax=1000,tune_rotation=True,p10=4.6e-4,max_iter=5):
    """
     p10 is the initial guess for distortion parameter, distortion = p10*rad^3
    """
    fac = 57.29577951

    base=fits_file.replace('.fits','')
    srcfile=base+'_dir/'+base+'_radec.txt.match.txt'

    dat = loadtxt(srcfile,usecols=(2,4,5,8,9)).T
    ii = dat[0].argsort()[:nmax]
    r0,d0,x,y = dat[1:,ii]

    hdr = getheader(fits_file)
    i,j = xy2ij(x,y,hdr,distort=False)
    ra0,dec0 = hdr['CRVAL1'],hdr['CRVAL2']
    c=cos(dec0/fac)

    sd0,cd0 = sin(d0/fac),cos(d0/fac)
    sr0,cr0 = sin(r0/fac),cos(r0/fac)
    def get_i0j0(i00,j00):
        """
          apply the tan projection to get i0,j0 from the reference r0,d0
        """
        da0, dd0 = (i00/c+ra0)/fac, (j00+dec0)/fac
        sda0,cda0 = sin(da0),cos(da0)
        sdd0,cdd0 = sin(dd0),cos(dd0)
        cfac = cr0*cda0+sr0*sda0
        cos_c = sdd0*sd0 + cdd0*cd0*cfac
        i0 = fac*cd0*(sr0*cda0-cr0*sda0) / cos_c
        j0 = fac*(cdd0*sd0-sdd0*cd0*cfac) / cos_c
        return i0,j0

    def get_p1(rad,rad0,p10=4.6e-4):
        """
          returns the radial dependence parameter p1, after sigma clipping
        """
        dis = rad0-rad
        rad3 = rad**3
        k = nsigma_clip( (dis-rad3*p10)/(1.+rad3*p10*1800.),5)*(rad0>0.5)
        ks = k.sum()
        if (ks>10):
            p1 = (dis[k]*rad3[k]).mean() / (rad3[k]**2).mean()
            p1 = ( p10 + p1*ks/100. )/(1. + ks/100.)
            if (p1<p10/2 or p1>2*p10): p1=p10
            chi2 = ( (dis[k]-rad3[k]*p1)**2 ).mean()
        else:
            p1 = p10
            chi2 = ( (dis-rad3*p1)**2 ).mean()
        return p1,chi2

    def fit_optcen(par,i,j):
        """
          optimize to determine the optical center (changing CRVAL1 CRVAL2)
        """
        i00,j00 = max(min(par[0]/360.,0.5),-0.5),max(min(par[1]/360.,0.5),-0.5)
        i0,j0 = get_i0j0(i00,j00)
        return get_p1( sqrt((i-i00)**2+(j-j00)**2), sqrt(i0**2+j0**2), p10=p10)[1]

    if (tune_rotation):
        h = i**2+j**2 < 2
        if (h.sum()<10): h = i**2+j**2 < 16.

    res = array([0.,0.])
    matr = array([[1.,0],[0,1]])
    for niter in xrange(max_iter):
        res = fmin(fit_optcen,res,args=(i,j),disp=False); res = fmin(fit_optcen,res,args=(i,j),disp=False)
        #
        i00,j00 = max(min(res[0]/360.,0.5),-0.5),max(min(res[1]/360.,0.5),-0.5)
        i0,j0 = get_i0j0(i00,j00)
        rad2 = (i-i00)**2+(j-j00)**2
        p1 = get_p1(sqrt(rad2),sqrt(i0**2+j0**2),p10=p10)[0]
        i1 = i + (i-i00)*rad2*p1
        j1 = j + (j-j00)*rad2*p1
        #
        if (tune_rotation):
            rii,rji = linfit(i1[h]-i00,i1[h]-i00-i0[h])[1], linfit(i1[h]-i00,j1[h]-j00-j0[h])[1]
            rjj,rij = linfit(j1[h]-j00,j1[h]-j00-j0[h])[1], linfit(j1[h]-j00,i1[h]-i00-i0[h])[1]
        else: rii=0;rij=0;rji=0;rjj=0;
        rsd1,rsd2 = abs(rji)*3600.,abs(rij)*3600.
        print """ Iter %d %.2e %.0f %.0f Resids: %.8f %.8f (arcsec) %s""" % (niter,p1,res[0]*10,res[1]*10,rsd1,rsd2,fits_file)
        #
        matr0 = array([[1-rii,-rij],[-rji,1-rjj]])
        i = i00 + matr0[0,0]*(i-i00) + matr0[0,1]*(j-j00)
        j = j00 + matr0[1,0]*(i-i00) + matr0[1,1]*(j-j00)
        matr = dot(matr0,matr)
        if (rsd1<1.e-3 and rsd2<1.e-3): break


    # diagnostic plot
    rad2 = (i-i00)**2+(j-j00)**2
    i1 = i + (i-i00)*rad2*p1
    j1 = j + (j-j00)*rad2*p1
    rad20 = i0**2+j0**2
    dis=sqrt(rad20)-sqrt(rad2)
    dism = sqrt(rad2)*rad2*p1
    dis1=sqrt(rad20)-sqrt((i1-i00)**2+(j1-j00)**2)

    plot (sqrt(rad2),dis*3600.,'bo',label='Uncorrected',alpha=0.5)
    s=rad2.argsort()
    plot (sqrt(rad2[s]),dism[s]*3600.,'k-',lw=3,alpha=0.5,label=("""%.1f""" % (p1*8*3600.))+r'$(r/2)^3$')
    plot (sqrt(rad2),dis1*3600.,'go',label='Corrected',alpha=0.5)

    k = nsigma_clip(dis1/(1.+dism*1800.),5)
    ylim((min(dis1[k].min(),dis[k].min())*3600,max(dis1[k].max(),dis[k].max())*3600.))
    grid()
    xlabel("Radius from Image Center [Degrees]",fontsize=16)
    ylabel("Source Offset [arcsec]",fontsize=16)
    legend(loc=2)
    savefig(fits_file.replace('.fits','_dist.jpg'))

    # p1*rad^3
    os.system("""sethead PV1_1=1.0 PV2_1=1.0 %s""" % fits_file)
    os.system("""sethead PV1_7=%.8f PV1_9=%.8f PV2_7=%.8f PV2_9=%.8f %s""" % (p1,p1,p1,p1,fits_file))
    hdr = getheader(fits_file)
    cd = array([[hdr['CD1_1'],hdr['CD1_2']],[hdr['CD2_1'],hdr['CD2_2']]])
    cd1 = dot(matr,cd)
    os.system("""sethead CD1_1=%.8f CD1_2=%.8f CD2_1=%.8f CD2_2=%.8f %s""" % (cd1[0,0],cd1[0,1],cd1[1,0],cd1[1,1],fits_file))
    hdr = getheader(fits_file)
    x,y = ij2xy(i00,j00,hdr,distort=True); x,y = round(x),round(y)
    r,d = xy2ad(x,y,hdr)
    os.system("""sethead CRPIX1=%.2f CRPIX2=%.2f CRVAL1=%.8f CRVAL2=%.8f %s""" % (x,y,r,d,fits_file))


if __name__ == "__main__":
    """
    """
    if (len(sys.argv)<2): usage()

    fits_file=sys.argv[1]
    if (os.path.exists(fits_file)==False): usage()

    p10=4.6e-4
    if (len(sys.argv)>2): p10=float(sys.argv[2])

    ddoti_distortion(fits_file,p10=p10)
