#!/usr/bin/python
"""
 sources2mask.py wt_file sources_file <outfile> <mask_fwhm>
  [requires mask.fits from sextractor]
"""
import os, sys
from pyfits import getdata,writeto
from scipy.ndimage import gaussian_filter
from scipy.signal import medfilt2d
from numpy import median,loadtxt,atleast_1d,sqrt

def usage():
    print __doc__
    sys.exit()


def sources2mask(wt_file,source_file,outfile='mask.fits',fwhm=7.):
    """
      Turn a sextractor mask.fits file into a 0,1 mask (0==masked),
    """

    print """Making object mask for: %s""" % source_file

    if (fwhm<1): fwhm = 1.
    sigma = fwhm/2.35

    # read in the weightmap
    wt=getdata(wt_file)
    wt_type = wt_file.split('.')[-2]
    if (wt_type=='rms'):
        dat_file = wt_file.replace('.rms.fits','.fits')
        rms = wt
    else:
        dat_file = wt_file.replace('.wt.fits','.fits')
        rms = 0.*wt
        h = wt>0
        rms[h] = 1./sqrt(wt[h])

    # make the source mask, from the sextractor source_file
    mask0,hdr = getdata(source_file,header=True)
    if (sigma>0):
        mask = gaussian_filter(mask0,sigma)
        mask[(mask0==0)*(mask<rms)]=0.
        mask[mask0>0]=1
        mask0 = 0.
    else:
        mask = mask0

    mask+=1; mask[mask!=1]=0

    if (sigma>0):
        mask = gaussian_filter(mask,sigma)
        mask[mask>0.9]=1
        mask[mask<=0.9]=0

    mask[wt==0]=0

    if (os.path.exists(outfile)): os.remove(outfile)
    writeto(outfile,mask,hdr)


if __name__ == '__main__':
    """
    """

    if (len(sys.argv)<3): usage()

    wtfile= sys.argv[1]
    if (os.path.exists(wtfile)==False): usage()
    sourcesfile= sys.argv[2]
    if (os.path.exists(sourcesfile)==False): usage()

    maskfile='mask.fits'
    if (len(sys.argv)>=3): maskfile=sys.argv[3]

    fwhm=3.0
    if (len(sys.argv)>4): fwhm=float(sys.argv[4])

    sources2mask(wtfile,sourcesfile,outfile=maskfile,fwhm=fwhm)
