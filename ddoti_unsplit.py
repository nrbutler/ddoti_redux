#!/usr/bin/python
"""
 ddoti_unsplit.py file_list [nsplit]
"""

import pyfits
from numpy import zeros,loadtxt
import sys,os

def usage():
    print __doc__
    sys.exit()

def ddoti_unsplit(file_list,nsplit=4):
    """
    split the input file into nsplit x nsplit
    """

    filename, file_extension = os.path.splitext(file_list)
    if (file_extension=='.fits'): files=[file_list]
    else: files = loadtxt(file_list,dtype='string').tolist()

    hdr0=pyfits.getheader('f00_'+files[0])
    a,b = hdr0['NAXIS2'],hdr0['NAXIS1']
    dat = zeros((a*nsplit,b*nsplit),dtype='float32')

    for file in files:

        hdr=pyfits.getheader('f00_'+file)
        for i in xrange(nsplit):
            for j in xrange(nsplit):
                file1="""f%d%d_%s""" % (i,j,file)
                dat[i*a:(i+1)*a,j*b:(j+1)*b] = pyfits.getdata(file1)

        pyfits.writeto(file,dat,hdr,clobber=True)


if __name__ == '__main__':
    """
    """
    if (len(sys.argv)<2): usage()

    file_list = sys.argv[1]

    nsplit=4
    if (len(sys.argv)>2): nsplit=int(sys.argv[2])

    ddoti_unsplit(file_list,nsplit=nsplit)
