#!/usr/bin/python
"""
 ddoti_unsplit.py fits_file [nsplit]
"""

from pyfits import getdata,getheader,writeto
from numpy import zeros
import sys,os

def usage():
    print __doc__
    sys.exit()

def ddoti_unsplit(fits_file,nsplit=4):
    """
    split the input file into nsplit x nsplit
    """

    hfile=fits_file.replace('.fits','.head')
    if (os.path.exists(hfile)):
        hdr0 = getheader(hfile)
        a0,b0 = hdr0['NAXIS2'],hdr0['NAXIS1']
        a,b = a0/nsplit,b0/nsplit
    else:
        hdr0 = getheader('f00_'+fits_file)
        a,b = hdr0['NAXIS2'],hdr0['NAXIS1']
        a0,b0 = nsplit*a,nsplit*b

    dat = zeros((a0,b0),dtype='float32')

    do_weight=False
    wfile0 = 'f00_'+fits_file.replace('.fits','.wt.fits')
    if (os.path.exists(wfile0)):
        wdat = zeros((a0,b0),dtype='float32')
        do_weight=True

    for i in xrange(nsplit):
        i1,i2 = i*a,(i+1)*a
        if (i==nsplit-1): i2=a0
        for j in xrange(nsplit):
            j1,j2 = j*b,(j+1)*b
            if (j==nsplit-1): j2=b0
            file1="""f%d%d_%s""" % (i,j,fits_file)
            dat[i1:i2,j1:j2] = getdata(file1)
            if (do_weight):
                wfile1=file1.replace('.fits','.wt.fits')
                wdat[i1:i2,j1:j2] = getdata(wfile1)

    writeto(fits_file,dat,hdr0,clobber=True)
    if (do_weight): writeto(fits_file.replace('.fits','.wt.fits'),wdat,hdr0,clobber=True)


if __name__ == '__main__':
    """
    """
    if (len(sys.argv)<2): usage()

    fits_file = sys.argv[1]

    nsplit=4
    if (len(sys.argv)>2): nsplit=int(sys.argv[2])

    ddoti_unsplit(fits_file,nsplit=nsplit)
