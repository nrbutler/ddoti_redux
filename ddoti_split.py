#!/usr/bin/python
"""
 ddoti_split.py file_list [nsplit]
"""

import pyfits
from numpy import median,zeros,loadtxt
import sys,os

def usage():
    print __doc__
    sys.exit()

def ddoti_split(file_list,nsplit=4,skip=10):
    """
    split the input file into nsplit x nsplit
    skip is for median calculation only
    """

    filename, file_extension = os.path.splitext(file_list)
    if (file_extension=='.fits'): files=[file_list]
    elif (file_extension=='.head'): files=[file_list]
    else: files = loadtxt(file_list,dtype='string',ndmin=1).tolist()

    hdr0=pyfits.getheader(files[0])
    a0,b0,c0,d0 = 1,hdr0['NAXIS1'],1,hdr0['NAXIS2']

    try: x0,y0 = hdr0['CRPIX1'],hdr0['CRPIX2']
    except: x0,y0 = b0/2.,d0/2.

    try: a,b,c,d = eval(hdr0['DATASEC'].replace('[','').replace(']','').replace(':',','))
    except: a,b,c,d = a0,b0,c0,d0

    sx0,sy0 = d-c+1,b-a+1
    sx,sy = sx0/nsplit,sy0/nsplit
    x = zeros((sx0,sy0),dtype='float32')

    do_bias=False
    try:
        a1,b1,c1,d1 = eval(hdr0['BIASSEC'].replace('[','').replace(']','').replace(':',','))
        do_bias=True
    except: pass

    for file in files:

        headfile=False
        filename, file_extension = os.path.splitext(file)
        if (file_extension=='.head'): headfile=True

        if (do_bias):
            hdu=pyfits.open(file); hdr=hdu[0].header
            try: bz=hdr['BZERO']
            except: bz=0.
            md0 = 0.
            if (b1>a1): md0 = median( hdu[0].section[c1-1:d1,a1-1:b1] )
            sky = median( hdu[0].section[d0/2-500:d0/2+500,b0/2-500:b0/2+500] ) - md0
            print """sethead SKYLEV=%.2f BZERO=%.2f %s""" % (sky,bz-md0,file)
            hdu.close()

        for i in xrange(nsplit):
            j0,j1 = i*sx+c,(i+1)*sx+c-1
            if (i==nsplit-1): j1 = sx0+c-1
            for j in xrange(nsplit):
                i0,i1 = j*sy+a,(j+1)*sy+a-1
                if (j==nsplit-1): i1 = sy0+a-1
                file1="""f%d%d_%s""" % (i,j,file)
                if (headfile): str1="""cp %s %s ; sethead I0=%d I1=%d J0=%d J1=%d NAXIS1=%d NAXIS2=%d CRPIX1=%s CRPIX2=%s %s\n""" % (file,file1,i0,i1,j0,j1,i1-i0+1,j1-j0+1,str(x0-j*sy-a+1),str(y0-i*sx-c+1),file1)
                else: str1="""getfits %s %d-%d %d-%d -o %s""" % (file,i0,i1,j0,j1,file1)
                print str1


if __name__ == '__main__':
    """
    """
    if (len(sys.argv)<2): usage()

    file_list = sys.argv[1]
    if (os.path.exists(file_list)==False): usage()

    nsplit=4
    if (len(sys.argv)>2): nsplit=int(sys.argv[2])

    ddoti_split(file_list,nsplit=nsplit)
