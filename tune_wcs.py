#!/usr/bin/python
"""
tune_wcs.py fitsfile calfile.txt
"""

import sys,os

from numpy import loadtxt,sqrt,array,dot,cos,sin,arange,vstack,median
from fit_wcs import ad2ij,xy2ij,xy2ad,ij2ad
from pyfits import getheader,getdata
from scipy.optimize import fmin
from nsigma_clip import nsigma_clip
from scipy.spatial import ckdtree

def usage():
     print __doc__
     sys.exit()


def tune_wcs(fits_file,catfile):
    """
    tune up a wcs solution in one file to match the catfile
    """

    hdr = getheader(fits_file)

    base=fits_file.replace('.fits','')
    radec_file=base+'_dir/'+base+'_radec.txt'
    if (os.path.exists(radec_file)==False):
        print """No radec file %s found!""" % radec_file
        sys.exit()

    phot_dat = loadtxt(radec_file,ndmin=2).T
    phot_dat = phot_dat[:,phot_dat[3]<999]

    ref_dat = getdata(catfile)

    idx0 = 1+arange(len(ref_dat[0]))
    h = ref_dat[3]<999
    ref_dat = ref_dat[:,h]
    idx0 = idx0[h]

    # match to the catalog using intermediate coordinates (i,j), using kdtree
    x,y = phot_dat[7:9]
    i,j = xy2ij(x,y,hdr,distort=True)

    ra0,dec0 = ref_dat[:2]
    i0,j0 = ad2ij(ra0,dec0,hdr)

    tre = ckdtree.cKDTree(vstack((i0,j0)).T)
    dis,ii = tre.query(vstack((i,j)).T)

    if (len(ii)>0):
        ref_dat = ref_dat[:,ii]
        i0,j0 = i0[ii],j0[ii]
        idx0 = idx0[ii]

        phot_dat[2] += median(ref_dat[2] - phot_dat[2])

        delta_i,delta_j = i-i0,j-j0
        k = nsigma_clip(delta_i,2) * nsigma_clip(delta_j,2)
        delta_i,delta_j = delta_i[k].mean(),delta_j[k].mean()
        i -= delta_i; j -= delta_j
        crval1,crval2 = ij2ad(-delta_i,-delta_j,hdr)

        par=array([0.,0.,0.])

        def myfun(par):
            a,b,th = par
            s,c = sin(th),cos(th)
            di = i0 - (1+a)*(c*i+s*j)
            dj = j0 - (1+b)*(c*j-s*i)
            k = nsigma_clip(di,5) * nsigma_clip(dj,5)
            return ( (di[k]*3600.)**2+(dj[k]*3600.)**2 ).mean()

        res = fmin(myfun,par,disp=0)
        res = fmin(myfun,res,disp=0)
        a,b,th = res
        s,c = sin(th),cos(th)

        matr = array([ [(1+a)*c,(1+a)*s],[-(1+b)*s,(1+b)*c] ] )
        c11,c12,c21,c22=hdr['CD1_1'],hdr['CD1_2'],hdr['CD2_1'],hdr['CD2_2']
        cd = array([ [c11,c12],[c21,c22] ])
        cd1 = dot(matr,cd)
        hdr['CD1_1'],hdr['CD1_2'],hdr['CD2_1'],hdr['CD2_2'],hdr['CRVAL1'],hdr['CRVAL2'] = cd1[0,0],cd1[0,1],cd1[1,0],cd1[1,1],crval1,crval2
        str="""sethead CD1_1=%.8f CD1_2=%.8f CD2_1=%.8f CD2_2=%.8f CRVAL1=%.8f CRVAL2=%.8f %s""" % (cd1[0,0],cd1[0,1],cd1[1,0],cd1[1,1],crval1,crval2,fits_file)
        os.system(str)

        # need to write out radec.txt.match.txt
        r,d = xy2ad(x,y,hdr,distort=True)

        matchfile=radec_file+'.match.txt'
        of=open(matchfile,'w')
        of.write("#RA Dec mag dmag RA_cat Dec_cat cat_mag cat_dmag x y xa ya x2a y2a expos num\n")
        for i in range(len(ref_dat[0])):
            ra,dec,mag0,dmag0 = ref_dat[:4,i]
            mag,dmag = phot_dat[2:4,i]
            x,y,xa,ya,x2a,y2a,expos,idx = phot_dat[7:,i]
            of.write("""%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n""" % (r[i],d[i],mag,dmag,ra,dec,mag0,dmag0,x,y,xa,ya,x2a,y2a,expos,idx0[i]))
        of.close()


if __name__ == "__main__":
    """
    """
    if (len(sys.argv)<3): usage()

    fits_file=sys.argv[1]
    if (os.path.exists(fits_file)==False): usage()

    catfile=sys.argv[2]
    if (os.path.exists(catfile)==False): usage()

    tune_wcs(fits_file,catfile)
