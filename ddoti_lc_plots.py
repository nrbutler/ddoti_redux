#!/usr/bin/python
"""
 ddoti_lc_plots.py datafile radecfile timefile t0
"""
import sys,os
from numpy import loadtxt,unique,floor,sqrt,median,ones,array,hstack,log10
from ut2gps import ut2gps

import matplotlib as mpl
mpl.use('Agg')

from matplotlib.pyplot import plot,xlabel,ylabel,savefig,ylim,xlim,errorbar,clf,scatter,title,legend

from linfit import linfit

def usage():
    print __doc__
    sys.exit()

def lc_plots_ddoti(datafile,radecfile,timefile,t0,nmin=2,nmax=100):
    """
       allow a maximum of nmax sources through
       if present not just in the stack, must be present at least nmin times
    """
    ofile=datafile+'.summary.txt'
    of = open(ofile,'w')
    of.write("#   id    RA           DEC          mag      dmag     fwhm     slope    dslope       chi2   n_detect\n")

    data = loadtxt(datafile,ndmin=2).T
    if (len(data[0])==0):
        of.close()
        sys.exit()

    id,mag,dmag,fwhm,expos,epoch,part = data
    id = id.astype('int32')
    epoch = epoch.astype('int32')

    gps0 = ut2gps(t0)
    times = loadtxt(timefile,dtype='string',usecols=(0,),ndmin=1)
    gps = ut2gps(times)-gps0

    expos = loadtxt(timefile,usecols=(1,),ndmin=1)
    gps0=0.
    if (gps.min()<expos.min()): gps0 += 3600.

    gps += expos/2. + gps0
    gps /= 3600.
    expos /= 3600.

    gps = hstack((gps.mean(),gps))
    expos = hstack((gps.max()-gps.min()+expos.mean(),expos))

    nepoch = len(unique(epoch))-1

    cat_dat = loadtxt(radecfile,ndmin=2).T
    if (len(cat_dat[0])==0):
        of.close()
        sys.exit()

    ra,dec = cat_dat[:2]
    remaining = ones(len(ra),dtype='bool')

    uid = unique(id)
    i0=1
    for id0 in uid[:nmax]:
        h = id==id0
        n,p,m,dm,f=epoch[h],part[h],mag[h],dmag[h],fwhm[h]

        # get rid of multiple detections in an epoch
        good = ones(len(n),dtype='bool')
        n1 = unique(n[n>0])
        for nn in n1:
            h = n==nn
            if (h.sum()>0): good[h] = m[h]<=m[h].min()

        if (good.sum()>1): n,p,m,dm,f = n[good],p[good],m[good],dm[good],f[good]

        h1 = (n>0)*(dm<999)
        n0 = h1.sum()
        mag0,dmag0 = m[n==0],dm[n==0]
        if (n0<nmin):
            if (nepoch<nmin):
                of.write("""%6d %12.8f %12.8f %8.4f %8.4f %8.4f %8s %8s %10s %4d\n""" % (i0,ra[id0-1],dec[id0-1],mag0,dmag0,f[n==0],'NA','NA','NA',n0))
                i0+=1
            else: remaining[id0-1] = False
            continue

        x,y = floor(p/10.),p%10
        d= sqrt( 0.5*((x/1.5-1)**2+(y/1.5-1)**2) )*100

        j=dm<999
        j1=j*(n>0)
        dm[~j] = 0.5; m[~j] += 0.5
        clf()
        errorbar(gps[n],m,xerr=expos[n]/2.,yerr=dm,marker='o',capsize=0,linestyle='None',markersize=0,mew=1)
        ylim((m.max()+0.5,m.min()-0.5))
        scatter(gps[n[j]],m[j],s=d[j])
        plot (gps[n[~j]],m[~j]+0.5,'bv')
        plot (gps[n[~j]],m[~j]-0.5,'bo')

        x=-2.5*log10(gps[n])
        res = linfit(x[j1],m[j1],dy=dm[j1])
        sig = (m[j1]-res[0]-res[1]*x[j1]).std()
        if (sig<0.01): sig=0.01
        var = dm[j1]**2 + sig**2

        m0 = (m[j1]/var).sum()/(1./var).sum()
        chi2 = ((m[j1]-m0)**2/var).mean()

        xm = x[j1].mean(); ym = m[j1].mean()
        res = linfit(x[j1]-xm,m[j1]-ym,dy=sqrt(var))
        slp=res[1]
        dslp=sqrt(res[2][1][1])
        # handle over-fitting
        if (j1.sum()<=2):
            if (slp<0): slp=min(slp+dslp,0.)
            else: slp=max(slp-dslp,0.)
        of.write("""%6d %12.8f %12.8f %8.4f %8.4f %8.4f %8.4f %8.4f %10.2f %4d\n""" % (i0,ra[id0-1],dec[id0-1],mag0,dmag0,f[n==0],slp,dslp,chi2,n0))
        ii = gps[n[j1]].argsort()
        plot (gps[n[j1]][ii],res[0]+slp*(x[j1][ii]-xm)+ym,label="""%.4f +/- %.4f""" % (slp,dslp))
        legend()

        dt0 = expos.min()
        xlim((gps.min()-dt0,gps.max()+dt0))
        xlabel("""Time Since %s - %.2f [hours]""" % (t0,gps0/3600.),fontsize=16)
        ylabel("USNO R",fontsize=16)
        title("""Light Curve for Source %d""" % i0)
        savefig("""lc_%d.jpg""" % i0)
        i0+=1

    of.close()

    cat_dat = cat_dat[:,remaining]

    of=open(radecfile,'r')
    ln1 = of.readline()
    of.close()
    of=open(radecfile,'w')
    if (ln1[0]=='#'): of.write(ln1)
    for i in xrange(len(cat_dat[0])):
        ra,dec,mag,dmag,mag_big,dmag_big,fwhm,x,y,xa,ya,x2a,y2a,expos,idx = cat_dat[:,i]
        of.write("""%f %f %f %f %f %f %f %f %f %f %f %f %f %f %d\n""" % (ra,dec,mag,dmag,mag_big,dmag_big,fwhm,x,y,xa,ya,x2a,y2a,expos,i+1))

    of.close()


if __name__ == "__main__":
    """
    make summary photometry plots
    """
    if (len(sys.argv)<5): usage()

    infile=sys.argv[1]
    if (os.path.exists(infile)==False): usage()
    radecfile=sys.argv[2]
    if (os.path.exists(radecfile)==False): usage()
    timefile=sys.argv[3]
    if (os.path.exists(timefile)==False): usage()
    t0=sys.argv[4]

    nmin=2
    if (len(sys.argv)>5): nmin=int(sys.argv[5])
    nmax=100
    if (len(sys.argv)>6): nmax=int(sys.argv[6])

    lc_plots_ddoti(infile,radecfile,timefile,t0=t0,nmin=nmin,nmax=nmax)
