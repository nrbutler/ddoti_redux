#!/usr/bin/python
"""
 weight2rms.py weight_image data_image rms_image
"""
import os, pyfits, sys
from numpy import sqrt,ones,zeros

def usage():
    print __doc__
    sys.exit()


def make_rms(data,weight,sky_levl=0.,wtmin=1.e-3):
    """
     make the rms from the data and weight map
    """
    rms = zeros(data.shape,dtype='float32')

    h = weight>wtmin
    rms[h] = sqrt( (1 + data[h].clip(0)/sky_levl)/weight[h] )

    return rms


def weight2rms(weight_file,data_file,rms_file):
    """
    Take a sextractor exposure weight file and transform to rms.
    """

    print """Making rms image: %s from weight image %s""" % (rms_file,weight_file)

    dat = pyfits.getdata(data_file)
    hdr = pyfits.getheader(data_file)
    wt = pyfits.getdata(weight_file)

    try: sky_levl=hdr['SKYLEV']
    except: sky_levl=0.

    rms = make_rms(dat,wt,sky_levl=sky_levl)

    if (os.path.exists(rms_file)): os.remove(rms_file)
    pyfits.writeto(rms_file,rms,hdr)


if __name__ == '__main__':
    """
    """

    if (len(sys.argv)<4): usage()

    weight_file= sys.argv[1]
    if (os.path.exists(weight_file)==False): usage()
    data_file= sys.argv[2]
    if (os.path.exists(data_file)==False): usage()
    rms_file= sys.argv[3]

    weight2rms(weight_file,data_file,rms_file)
