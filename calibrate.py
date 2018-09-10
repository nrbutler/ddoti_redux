#!/usr/bin/python
"""
calibrate.py radec.txt catalog_radec.txt <phot_outfile.txt> <cat_outfile.txt>
   for comparing a sextractor radec list to catalogs
"""
import sys,os
from quick_mode import quick_mode
from robust_mean import robust_mean
from numpy import loadtxt,sqrt,log10,arange,where,median,zeros

def usage():
    print __doc__
    sys.exit()


def calibrate(infile='radec.txt',calfile='cat_radec.txt',outfile='photometry.txt',matchfile='match.txt',sys_err=0.001,cal_mag_max=16.0,ab_offset=0.58):

    phot_dat=loadtxt(infile,ndmin=2).T
    idx = phot_dat[-1].astype('int32')
    gphot = (phot_dat[2]<999)*(phot_dat[3]<999)*(idx>0)
    idx = idx[gphot]
    mag,dmag = phot_dat[2:4,gphot]; mag_big,dmag_big = phot_dat[4:6,gphot]; fwhm = phot_dat[6,gphot]

    lr1='#\n'
    f=open(infile,'r')
    try:
        lr1=f.readline(); dat=lr1.split();
        gain, dt, am, sex_zero, t0, t1 = float(dat[-12]), float(dat[-10]), float(dat[-8]), float(dat[-6]), dat[-4], dat[-2]
    except:
        gain,dt,am=1.,1.,1.
        sex_zero = 25.0
        t0, t1 = "0.0","0.0"
    f.close()

    # insert x,y dependent systematic error
    dmag = sqrt( dmag**2 + sys_err**2 )

    cdat=loadtxt(calfile,usecols=(2,3,4),ndmin=2).T
    cdat = cdat[:,(cdat[1]<999)*(cdat[2]<999)*(cdat[1]<cal_mag_max)]
    idx0 = cdat[0].astype('int32')

    # w matchlist refers to index in full star list
    # idx and idx0 don't necessarily have all those
    ii0 = zeros(max(idx0.max(),idx.max()),dtype='int32')-1
    ii0[idx0-1] = arange(len(idx0))
    w = ii0[idx-1]

    w0=w>=0
    nn = w0.sum()
    if (nn>0):

        print (" Matches %d" % nn)
        print (""" Time Range %s - %s""" % (t0,t1))

        w = w[w0]
        idx0 = idx0[w]
        mag0,dmag0 = cdat[1:,w]

        ii = where(gphot)[0]
        gphot[ii[~w0]] = False

        mag,dmag,mag_big,dmag_big,fwhm = mag[w0],dmag[w0],mag_big[w0],dmag_big[w0],fwhm[w0]

        fwhm0 = quick_mode(fwhm)
        print (""" Median Sextractor FWHM %.2f""" % fwhm0)

        mag_offset = mag0-mag+ab_offset
        dmag_offset = sqrt(dmag0**2+dmag**2)
        mag_offset, dmag_offset = robust_mean(mag_offset,dmag_offset)
        print (""" Magnitude Offset %.4f +/- %0.4f""" % (mag_offset,dmag_offset))

        phot_dat[2] += mag_offset
        phot_dat[4] += mag_offset
        h = phot_dat[3]<999
        phot_dat[3,h] = sqrt( phot_dat[3,h]**2 + dmag_offset**2 + sys_err**2 )

        # calculate the aperture correction
        ap_corr = median(mag_big-mag)
        zero_pt = sex_zero + mag_offset + 2.5*log10(gain/dt)+am - ap_corr
        print (""" Median Zero Point %.3f [gain=%.2f, dt=%.1f, am_corr=%.3f]""" % (zero_pt,gain,dt,am))
        print (""" Sextractor Aperture Correction (mag) %.3f""" % ap_corr)

        mag10s = quick_mode( mag-2.5*log10(dmag/0.10857) )
        print (" 10-sigma limiting magnitude %.2f" % (mag_offset+mag10s))

        of = open(outfile,'w')
        alr1 = lr1.split()
        alr1[-6] = """%f""" % (sex_zero+mag_offset)
        of.write(' '.join(alr1)+'\n')
        for i in xrange(len(phot_dat[0])):
            (ra,dec,mag,dmag,mag_big,dmag_big,fwhm,x,y,xa,ya,x2a,y2a,expos,idx) = phot_dat[:,i]
            of.write("""%f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n""" % (ra,dec,mag,dmag,mag_big,dmag_big,fwhm,x,y,xa,ya,x2a,y2a,expos,idx))

        of.close()

        f=open(matchfile,'w')

        f.write("#RA Dec mag dmag RA_cat Dec_cat cat_mag cat_dmag x y xa ya x2a y2a expos num\n")
        f.write("## Time Range %s - %s\n""" % (t0,t1))
        f.write("## Median Sextractor FWHM %.2f\n" % fwhm0)
        f.write("## Sextractor Mag Zero Point = %.4f +/- %.4f\n" % (sex_zero+mag_offset,dmag_offset))
        f.write("## Median Zero Point %.2f\n" % zero_pt)
        f.write("## 10-sigma limiting magnitude %.2f\n" % (mag_offset+mag10s))
        f.write("## Magnitude Offset %.4f +/- %0.4f\n" % (mag_offset,dmag_offset))

        phot_dat = phot_dat[:,gphot]
        for j in range(len(mag0)):
            (ra,dec,mag,dmag,mag_big,dmag_big,fwhm,x,y,xa,ya,x2a,y2a,expos,idx) = phot_dat[:,j]
            f.write("""%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n""" % (ra,dec,mag,dmag,ra,dec,mag0[j],dmag0[j],x,y,xa,ya,x2a,y2a,expos,idx0[j]))

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

    cal_mag_max=16.0
    if (len(sys.argv)>5): cal_mag_max=float(sys.argv[5])

    ab_offset=0.58
    if (len(sys.argv)>6): ab_offset=float(sys.argv[6])

    calibrate(infile,calfile,outfile,matchfile,cal_mag_max=cal_mag_max,ab_offset=ab_offset)
