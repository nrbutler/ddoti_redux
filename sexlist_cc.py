#!/usr/bin/python
"""
sexlist_cc.py fits_file shift_cat_file ref_cat_file cam
"""

import sys,os
from numpy import abs,where,sqrt,loadtxt,log,zeros,array,hstack,std,newaxis,median,pi,cos
from nsigma_clip import nsigma_clip
try: from weave import inline
except: from scipy.weave import inline
from coord import match
import pyfits
from robust_mean import robust_mean
from fit_wcs import ad2ij,ij2ad

def usage():
    print __doc__
    sys.exit()


def sexlist_cc(fits_file,shift_file,ref_file,x0=0.,y0=0.,search_rad=0.,cam='X',mag_diff_max=25.,score_min=-4,sigma=2.,nmatch_min=2,rad_max=0.):
    """
      Finds the offset in x and y between two images.  x and y are
      the intermediate WCS coordinates, which take into account 
      rotations, pixel scale, and distortions.

      assumes 2 images of the same WCS orientation, only differing by
      the cetral RA, Dec (CRVAL1,CRVAL2)

      search_rad: pixel (defaults to 0.75 imsize)
      sigma: psf sigma in arcsec (default 2)
      x0,y0: initial pixel offset guess
      rad_max: only used image sources within this radius from center

      Prints lists of image offsets.
    """
    hdr=pyfits.getheader(fits_file)

    try: platescale=sqrt( hdr['CD1_1']**2+hdr['CD1_2']**2 )*3600.
    except: platescale=0.3

    if (search_rad==0):
        try: search_rad = 0.75*max(hdr['NAXIS1'],hdr['NAXIS2'])
        except: search_rad = 100.

    # we work in pixels
    sigma /= platescale
    if (sigma<0.5): sigma=0.5

    nextend=0
    try:
        nextend = hdr['NEXTEND']
        hdr=pyfits.getheader(fits_file,1)
    except: pass

    using_stack_ref=False
    if (ref_file[0:6]=='stack_'): using_stack_ref=True

    using_stack_shf=False
    if (shift_file[0:6]=='stack_'): using_stack_shf=True

    newlist=False
    if (os.path.exists(ref_file+'.new.txt') and using_stack_ref==False):
        (x1,y1,m1,dm1,score) = loadtxt(ref_file+'.new.txt',unpack=True,usecols=(0,1,2,3,4),ndmin=2)
        newlist=True
        score = score.astype('int8')
        ii = score>=score_min
        m1 = m1[ii]; x1 = x1[ii]; y1 = y1[ii]; dm1 = dm1[ii]; score = score[ii]
    else:
        ref_dat = loadtxt(ref_file,ndmin=2).T
        (x1,y1,m1,dm1) = 1.*ref_dat[:4]
        score = zeros(len(m1),dtype='int8')
        x1,y1 = ad2ij(x1,y1,hdr)
        x1 *= 3600./platescale
        y1 *= 3600./platescale

    if os.path.exists(shift_file+'.new.txt'):
        (x2,y2,m2,dm2,score2) = loadtxt(shift_file+'.new.txt',unpack=True,usecols=(0,1,2,3,4),ndmin=2)
        ii = score2>=score_min
        m2 = m2[ii]; x2 = x2[ii]; y2 = y2[ii]; dm2 = dm2[ii]
    else:
        phot_dat = loadtxt(shift_file,ndmin=2).T
        (x2,y2,m2,dm2) = 1.*phot_dat[:4]
        x2,y2 = ad2ij(x2,y2,hdr)
        x2 *= 3600./platescale
        y2 *= 3600./platescale
        if (rad_max>0):
            h = x2**2+y2**2<rad_max**2
            x2,y2,m2,dm2 = x2[h],y2[h],m2[h],dm2[h]
        ii = x2.argsort()
        x2,y2,m2,dm2 = x2[ii],y2[ii],m2[ii],dm2[ii]

    h = (x1>x2.min()-search_rad)*(x1<x2.max()+search_rad)*(y1>y2.min()-search_rad)*(y1<y2.max()+search_rad)
    if (h.sum()>0):
        x1,y1,m1,dm1 = x1[h],y1[h],m1[h],dm1[h]

    n1,n2 = len(x1),len(x2)
    x2 -= x0; y2 -= y0

    # renormalize the fluxes
    f1 = -m1; f2 = -m2
    fm = f1.min(); f2m = f2.min()
    if (f2m<fm): fm=f2m
    f1=f1-fm+0.1; f2=f2-fm+0.1
    f1 /= median(f1); f2 /= median(f2)

    search_code = """
    unsigned long i,j,j0;
    int dx,dy;
    double fac;
    for (i=0;i<n1;i++) {
        for (j0=0;j0<n2;j0++) {
            if( (int)(x2a[j0]-x1a[i]+0.5)>=0 ) break;
        }
        for (j=j0;j<n2;j++) {
            dx = (int)(x2a[j]-x1a[i]+0.5);
            if (dx>=n) break;
            dy = (int)(y2a[j]-y1a[i]+0.5);
            if (dy>=0 && dy<n) {
                fac = (f1[i]+f2[j])/2.;
                cc[dy+n*dx] += f1[i]*f2[j]/(fac*fac);
            }
        }
    }
    """
    #
    # find matches at a resolution of 2/fwhm/sqrt(2)
    zoom_fac=0.6/sigma
    isr = int(round(search_rad*zoom_fac))
    x1a,y1a,x2a,y2a = x1*zoom_fac-isr,y1*zoom_fac-isr,x2*zoom_fac,y2*zoom_fac

    n=2*isr+1
    cc = zeros((n,n),dtype='float32')
    inline(search_code,['cc','x1a','y1a','x2a','y2a','f1','f2','n1','n2','n'])

    x2 += x0; y2 += y0

    k=where(cc==cc.max())
    delta_x, delta_y = (k[0][0]-isr)/zoom_fac+x0, (k[1][0]-isr)/zoom_fac+y0

    f= platescale/3600.
    x2p,y2p = x2-delta_x,y2-delta_y
    ii=match(x1*f,y1*f,x2p*f,y2p*f,3*sigma*platescale)
    w=where(ii['sep']>=0)

    nmatch = len(w[0])
    matched=False
    delta_mag = 0.

    if (nmatch>=nmatch_min):
        matched=True

        w1 = ii['ind'][w]
        delta_x, delta_y = x2[w1]-x1[w], y2[w1]-y1[w]
        k = nsigma_clip(delta_x,2) * nsigma_clip(delta_y,2)
        delta_x, delta_y = delta_x[k].mean(),delta_y[k].mean()

        delta_mag, ddelta_mag = robust_mean(m1[w][k]-m2[w1][k],sqrt(dm1[w][k]**2+dm2[w1][k]**2+0.003))
        if (abs(delta_mag)>mag_diff_max and using_stack_ref==False):
            sys.stderr.write("""Out-of-bounds magnitude difference %f for %s %s, bad cross-match\n""" % (delta_mag,shift_file,ref_file))
            matched=False
    else:
        sys.stderr.write("""Not enough matches found (%d<%d)\n""" % (nmatch,nmatch_min))

    if (matched):

        if (using_stack_shf==False and cam!='X'):

            # track matched and unmatched sources which could have been matched
            x1min,x1max = x1.min(),x1.max()
            y1min,y1max = y1.min(),y1.max()
            x2min,x2max = x2.min(),x2.max()
            y2min,y2max = y2.min(),y2.max()
            x2p = x2-delta_x; y2p = y2-delta_y
            x1p = x1+delta_x; y1p = y1+delta_y

            shf_in_ref = (x2p>=x1min)*(x2p<=x1max)*(y2p>=y1min)*(y2p<=y1max)
            ref_in_shf = (x1p>=x2min)*(x1p<=x2max)*(y1p>=y2min)*(y1p<=y2max)

            score_ref = zeros(len(x1),dtype='int8')
            score_ref[ ref_in_shf ] = -1; score_ref[w] = 1

            f=open(ref_file+'.new.txt','w')
            for i in xrange(len(x1)):
                str="""%f %f %f %f %d\n"""
                f.write(str % (x1[i],y1[i],m1[i],dm1[i],score_ref[i]+score[i]))

            # new sources not in the reference? save them for addition to the reference
            h = where( ~shf_in_ref )
            if (len(h[0])>0):
                for i in xrange(len(h[0])):
                    ii0=h[0][i]
                    dis2 = (x2p[ii0]-x1)**2 + (y2p[ii0]-y1)**2
                    if (dis2.min() > 10*sigma**2):
                        str="""%f %f %f %f 0\n"""
                        f.write(str % (x2p[ii0],y2p[ii0],m2[ii0],dm2[ii0]))

            f.close()

        dx,dy = -delta_x*platescale/3600., -delta_y*platescale/3600.
        ra1,dec1 = ij2ad(dx,dy,hdr)

        # correct the photometry file
        r,d = phot_dat[:2]
        i,j = ad2ij(r,d,hdr)
        if (hdr['CTYPE1'] == 'DEC--TAN'): hdr['CRVAL2'],hdr['CRVAL1'] = ra1,dec1
        else: hdr['CRVAL1'],hdr['CRVAL2'] = ra1,dec1
        r[:],d[:] = ij2ad(i,j,hdr)

        of=open(shift_file,'r')
        ln1 = of.readline()
        of.close()
        of=open(shift_file,'w')
        if (ln1[0]=='#'): of.write(ln1)
        for i in xrange(len(phot_dat[0])):
            ra,dec,mag,dmag,mag_big,dmag_big,fwhm,x,y,xa,ya,x2a,y2a,expos,idx = phot_dat[:,i]
            of.write("""%f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n""" % (ra,dec,mag,dmag,mag_big,dmag_big,fwhm,x,y,xa,ya,x2a,y2a,expos,idx))
        of.close()

        i,j = ad2ij(r,d,hdr)
        ra0,dec0,mag0,dmag0 = ref_dat[:4]
        i0,j0 = ad2ij(ra0,dec0,hdr)

        ii=match(i0,j0,i,j,10*sigma*platescale)
        w=where(ii['sep']>=0)[0]

        nn = len(w)
        if (nn>0):
            w1 = ii['ind'][w]
            idx0 = ref_dat[-1]
            ra0,dec0,mag0,dmag0,idx0 = ra0[w],dec0[w],mag0[w],dmag0[w],idx0[w]
            of=open(shift_file+'.match.txt','w')
            of.write("#RA Dec mag dmag RA_cat Dec_cat cat_mag cat_dmag x y xa ya x2a y2a expos num\n")
            for i in range(nn):
                j=w1[i]
                ra,dec,mag,dmag,mag_big,dmag_big,fwhm,x,y,xa,ya,x2a,y2a,expos,idx = phot_dat[:,j]
                of.write("""%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n""" % (ra,dec,mag,dmag,ra0[i],dec0[i],mag0[i],dmag0[i],x,y,xa,ya,x2a,y2a,expos,idx0[i]))
            of.close()

        if ( using_stack_shf==False or cam=='X' ):
            if (hdr['CTYPE1'] == 'DEC--TAN'):
                os.system("""sethead CRVAL2=%.8f CRVAL1=%.8f %s""" % (ra1,dec1,fits_file))
                os.system("""sethead CRVAL2=%.8f CRVAL1=%.8f %s -x ' ' 2>/dev/null""" % (ra1,dec1,fits_file))
            else:
                os.system("""sethead CRVAL1=%.8f CRVAL2=%.8f %s""" % (ra1,dec1,fits_file))
                os.system("""sethead CRVAL1=%.8f CRVAL2=%.8f %s -x ' ' 2>/dev/null""" % (ra1,dec1,fits_file))

    return delta_x,delta_y,nmatch,matched


def main():

    if (len(sys.argv)<4): usage()

    fits_file=sys.argv[1]
    if (os.path.exists(fits_file)==0): usage()
    shift_file=sys.argv[2]
    if (os.path.exists(shift_file)==0): usage()
    ref_file=sys.argv[3]
    if (os.path.exists(ref_file)==0): usage()
    cam = 'X'
    if (len(sys.argv)>4): cam = sys.argv[4]
    search_rad=0.
    if (len(sys.argv)>5): search_rad=float(sys.argv[5])
    x0=0.
    if (len(sys.argv)>6): x0=float(sys.argv[6])
    y0=0.
    if (len(sys.argv)>7): y0=float(sys.argv[7])
    mag_diff_max = 25.
    if (len(sys.argv)>8): mag_diff_max=float(sys.argv[8])
    nmatch_min=2
    if (len(sys.argv)>9): nmatch_min=int(sys.argv[9])
    rad_max=0.
    if (len(sys.argv)>10): rad_max=float(sys.argv[10])

    if (abs(x0)>1.e4 or abs(y0)>1.e4): x0,y0=0.,0.

    dx,dy,nmatch,state = sexlist_cc(fits_file,shift_file,ref_file,search_rad=search_rad,cam=cam,x0=x0,y0=y0,mag_diff_max=mag_diff_max,nmatch_min=nmatch_min,rad_max=rad_max)
    if (state):
        print "%.3f %.3f %.3f %.3f %d %s" % (x0,y0,dx,dy,nmatch,shift_file)
    else:
        sys.stderr.write("""Cross-match not found for %s %s, x0,y0= %f,%f (best_fit %d matches at %f %f)\n""" % (shift_file,ref_file,x0,y0,nmatch,dx,dy))


if __name__ == "__main__":
    main()
