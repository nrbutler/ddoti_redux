#!/usr/bin/python
"""
calibrate.py radec.txt catalog_radec.txt <phot_outfile.txt> <cat_outfile.txt>
   for comparing a sextractor radec list to catalogs
"""
import sys,os
from quick_mode import quick_mode
from robust_mean import robust_mean
from numpy import loadtxt,sqrt,where,log10,arange
from numpy import isnan,array,median,cos,pi

def usage():
    print __doc__
    sys.exit()


def calibrate(infile='radec.txt',calfile='cat_radec.txt',outfile='photometry.txt',matchfile='match.txt',sys_err=0.001):

    phot_dat=loadtxt(infile).T
    j = phot_dat[3]<99
    (ra,dec,mag,dmag,mag_big,dmag_big,fwhm,x,y,xa,ya,x2a,y2a,expos,idx)=phot_dat[:,j]
    idx = idx.astype('int16')

    lr1='#\n'
    f=open(infile,'r')
    try:
        lr1=f.readline(); dat=lr1.split();
        gain, dt, am, sex_zero, t0, t1 = float(dat[-12]), float(dat[-10]), float(dat[-8]), float(dat[-6]), dat[-4], dat[-2]
    except:
        gain,dt,am=1.,1.,1.
        sex_zero = 25.0
        t0, t1 = "0.0","0.0"

    lr2='#\n'
    maglim,ap_corr=999,999
    try:
        lr2=f.readline(); dat=lr2.split();
        if (dat[0]=='#'):
            maglim=float(dat[4])
            ap_corr=float(dat[-3])
    except:
        maglim=999
        ap_corr=999

    f.close()

    # insert x,y dependent systematic error
    dmag = sqrt( dmag**2 + sys_err**2 )
    phot_err = sqrt( dmag**2 + 0.01*(xa**2+ya**2) )

    (ra0,dec0,mag0,dmag0)=loadtxt(calfile,unpack=True,usecols=(0,1,2,3))
    idx0 = arange(len(ra0))
    j = dmag0<999
    ra0=ra0[j]; dec0=dec0[j]; mag0=mag0[j]; dmag0=dmag0[j]; idx0 = idx0[j]

    w0 = idx>0
    nn = w0.sum()
    w = idx[w0]-1
    if (nn>0):

        mag0 = mag0[w]; dmag0 = dmag0[w];
        ra0 = ra0[w]; dec0=dec0[w]
        idx0 = idx0[w]

        # calculate the aperture correction
        if (ap_corr>990): ap_corr,dap_corr = robust_mean(mag_big-mag,phot_err)

        dRA = median((ra[w0]-ra0)*cos(dec0[0]*pi/180.)*3600.)
        dDEC = median((dec[w0]-dec0)*3600.)

        fwhm0 = quick_mode(fwhm[w0])

        print (" Matches %d" % nn)
        print (""" Time Range %s - %s""" % (t0,t1))
        print (""" Median Sextractor FWHM %.2f""" % fwhm0)
        print (""" Position Offset (arcsec): %.2f, %.2f (%.2f)""") % (dRA,dDEC,sqrt(dRA**2+dDEC**2))

        print """ Correction Terms: 0 0 0 0"""

        mag_offset = mag0-mag[w0]
        dmag_offset = sqrt(dmag0**2+dmag[w0]**2)
        mag_offset, dmag_offset = robust_mean(mag_offset,dmag_offset)
        print (""" Magnitude Offset %.4f +/- %0.4f""" % (mag_offset,dmag_offset))
        dmag[dmag<=0] = 1.e-6

        # last term is the aperture correction from sextractor
        zero_pt = sex_zero + mag_offset + 2.5*log10(gain/dt)+am - ap_corr

        print (""" Median Zero Point %.3f [gain=%.2f, dt=%.1f, am_corr=%.3f]""" % (zero_pt,gain,dt,am))
        print (""" Sextractor Aperture Correction (mag) %.3f""" % ap_corr)

        mag += mag_offset
        mag_big += mag_offset
        h = dmag<999
        dmag[h] = sqrt( dmag[h]**2+dmag_offset**2 )

        diff=mag[w0]-mag0
        diff_err = diff/sqrt(dmag[w0]**2+dmag0**2)
        m0 = median(diff)
        m0e = median(diff_err)
        dm0 = 1.48*median(abs(diff-m0))
        dm0e = 1.48*median(abs(diff_err-m0e))
        j = abs(diff_err-m0e)<=2*dm0e
        diff = diff[j]
        diff_err = diff_err[j]
        m0 = median(diff)
        m0e = median(diff_err)
        dm0 = 1.48*median(abs(diff-m0))
        dm0e = 1.48*median(abs(diff_err-m0e))

        if (isnan(dm0) or isnan(dm0e)): dm0,dm0e=0.,0.

        print (""" Median scatter %.3f (%.2f)""" % (dm0,dm0e))
        j=where( abs(mag[w0]-mag0)<sqrt(dmag[w0]**2+dmag0**2) )
        j1=where( abs(mag[w0]-mag0)<3.*sqrt(dmag[w0]**2+dmag0**2) )
        print (""" Fraction within 1-sigma/3-sigma: %.2f""" % (1.*len(j[0])/len(j1[0])))

        mag10s = quick_mode( mag[w0]-2.5*log10(dmag[w0]/0.10857) )
        print (" 10-sigma limiting magnitude %.2f" % mag10s)
        print (" Maglim 10-sigma %.4f - 2.5*log10(rms)" % (mag_offset+maglim))

        of = open(outfile,'w')
        alr1 = lr1.split()
        alr1[-6] = """%f""" % (sex_zero+mag_offset)
        of.write(' '.join(alr1)+'\n'); of.write(lr2)
        for i in xrange(len(mag)):
            of.write("""%f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n""" % (ra[i],dec[i],mag[i],dmag[i],mag_big[i],dmag_big[i],fwhm[i],x[i],y[i],xa[i],ya[i],x2a[i],y2a[i],expos[i],idx[i]))

        of.close()

        f=open(matchfile,'w')

        f.write("#RA Dec mag dmag RA_cat Dec_cat cat_mag cat_dmag x y xa ya x2a y2a expos num\n")
        f.write("## Time Range %s - %s\n""" % (t0,t1))
        f.write("## Median Sextractor FWHM %.2f\n" % fwhm0)
        f.write("## Sextractor Mag Zero Point = %.4f +/- %.4f\n" % (sex_zero+mag_offset,dmag_offset))
        f.write("## Correction Terms: 0 0 0 0\n")
        f.write("## Median Zero Point %.2f\n" % zero_pt)
        f.write("## 10-sigma limiting magnitude %.2f\n" % mag10s)
        f.write("## Magnitude Offset %.4f +/- %0.4f\n" % (mag_offset,dmag_offset))

        ra,dec,mag,dmag,x,y,xa,ya,x2a,y2a,expos = ra[w0],dec[w0],mag[w0],dmag[w0],x[w0],y[w0],xa[w0],ya[w0],x2a[w0],y2a[w0],expos[w0]
        for j in range(len(mag0)):
            f.write("""%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n""" % (ra[j],dec[j],mag[j],dmag[j],ra0[j],dec0[j],mag0[j],dmag0[j],x[j],y[j],xa[j],ya[j],x2a[j],y2a[j],expos[j],idx0[j]))

        f.close()

    else:
        print "No matches!"


if __name__ == "__main__":

    if (len(sys.argv)<3): usage()

    infile=sys.argv[1]
    if (os.path.exists(infile)==0): usage()
    calfile=sys.argv[2]
    if (os.path.exists(calfile)==0): usage()

    outfile=infile+'.photometry.txt'
    if (len(sys.argv)>3): outfile=sys.argv[3]

    matchfile=infile+'.match.txt'
    if (len(sys.argv)>4): matchfile=sys.argv[4]

    calibrate(infile,calfile,outfile,matchfile)
