"""Astronomical coordinate functions."""
import numpy as np
from numpy import arccos,sin,cos
from math import pi
# constants
DEG_PER_HR = 360. / 24.             # degrees per hour
DEG_PER_MIN = DEG_PER_HR / 60.      # degrees per min
DEG_PER_S = DEG_PER_MIN / 60.       # degrees per sec
DEG_PER_AMIN = 1./60.               # degrees per arcmin
DEG_PER_ASEC = DEG_PER_AMIN / 60.   # degrees per arcsec
RAD_PER_DEG = pi / 180.             # radians per degree

def ang_sep(ra1, dec1, ra2, dec2):
    """ Returns the angular separation (in units of degrees) on the
    celestial sphere between two ra/dec coordinates given in degrees.
    Accepts numpy arrays.

    Note: only works for separations larger than about 0.1 arcsec.
    Smaller separations are always returned as zero, because of
    floating point effects.

    >>> np.allclose(ang_sep(2, 0, 4, 0), 2.)
    True
    >>> np.allclose(ang_sep(359, 0, 1, 0), 2.)
    True
    >>> np.allclose(ang_sep(0, 20, 0, -10), 30)
    True
    >>> np.allclose(ang_sep(7, 20, 8, 40), 20.018358)
    True
    >>> np.allclose(ang_sep(7, 20, 250, -50.3), 122.388401)
    True
    >>> ang_sep(-1, 10, 240, -10)           # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    ValueError: RA outside sensible limits. -1
    >>> ras = [24.5,23.6]; decs = [66.89,67.01]
    >>> ref_sep = np.array([3.520032,  3.26675])
    >>> np.allclose(ang_sep(20. ,70., ras, decs), ref_sep)
    True
    >>> ras = [24.5, 23.6]; decs = [91.89, 67.01]
    >>> ang_sep(20.0, 70.0, ras, decs)      # doctest: +ELLIPSIS
    Traceback (most recent call last):
    ...
    ValueError: Dec outside sensible limits. [ 91.89  67.01]
    """ 
    ra1 = np.asarray(ra1);  ra2 = np.asarray(ra2)
    dec1 = np.asarray(dec1);  dec2 = np.asarray(dec2)
    # error checking
    # for ra in ra1,ra2:
    #     if not ((0. <= ra)&(ra < 360.)).all():
    #         raise ValueError("RA outside sensible limits. %s" % ra)
    # for dec in dec1,dec2:
    #     if np.any(np.abs(dec) > 90.):
    #         raise ValueError("Dec outside sensible limits. %s" % dec)

    ra1 = ra1 * RAD_PER_DEG           # convert to radians
    ra2 = ra2 * RAD_PER_DEG
    dec1 = dec1 * RAD_PER_DEG
    dec2 = dec2 * RAD_PER_DEG

    sra1 = sin(ra1);  sra2 = sin(ra2)
    cra1 = cos(ra1);  cra2 = cos(ra2)
    sdec1 = sin(dec1);  sdec2 = sin(dec2)
    cdec1 = cos(dec1);  cdec2 = cos(dec2)

    csep = cdec1*cdec2*(cra1*cra2 + sra1*sra2) + sdec1*sdec2

    # An ugly work-around for floating point issues.
    #if np.any(csep > 1):  print csep
    csep = np.where(csep > 1., 1., csep)

    degsep = arccos(csep) / RAD_PER_DEG
    # only works for separations > 0.1 of an arcsec or  >~2.7e-5 dec
    degsep = np.where(degsep < 1e-5, 0, degsep)
    return degsep


def match(ra1, dec1, ra2, dec2, tol, allmatches=False):
    """
    match(ra1, dec1, ra2, dec2, tol)

    Given two sets of numpy arrays of ra,dec and a tolerance tol
    (float), returns an array of indices and separations with the same
    length as the first input array.  If and index is > 0, it is the
    index of the closest matching second array element within tol
    arcsec.  If it's -1, then there was no matching ra/dec within tol
    arcsec.

    if allmatches = True, then for each object in the first array,
    return the index and separation of everything in the second array
    within the search tolerance, not just the closest match.

    Note to get the indices of objects in ra2, dec2 without a match, use

    imatch = match(ra1, dec1, ra2, dec2, 2.)
    inomatch = numpy.setdiff1d(np.arange(len(ra2)), set(imatch))

    """
    from numpy.core.records import fromarrays
    
    ra1,ra2,dec1,dec2 = map(np.asarray, (ra1, ra2, dec1, dec2))

    abs = np.abs

    isorted = ra2.argsort()
    sdec2 = dec2[isorted]
    sra2 = ra2[isorted]

    LIM = tol * DEG_PER_ASEC

    match = []
    # use mean dec, assumes decs similar
    decav = np.mean(sdec2.mean() + dec1.mean())
    RA_LIM = LIM / cos(decav * RAD_PER_DEG)

    for ra,dec in zip(ra1,dec1):
        i1 = sra2.searchsorted(ra - RA_LIM)
        i2 = i1 + sra2[i1:].searchsorted(ra + RA_LIM)
        #print i1,i2
        close = []
        for j in xrange(i1,i2):
            if abs(dec - sdec2[j]) > LIM:
                continue
            else:
                # if ras and decs are within LIM arcsec, then
                # calculate actual separation:
                disq = ang_sep(ra, dec, sra2[j], sdec2[j])
                close.append((disq, j))

        close.sort()
        if not allmatches:
            # Choose the object with the closest separation inside the
            # requested tolerance, if one was found.
            if len(close) > 0:
                min_dist, jmin = close[0]
                if min_dist < LIM:
                    match.append((isorted[jmin], min_dist))
                    continue
            # otherwise no match
            match.append((-1,-1.))
        else:
            # append all the matching objects
            jclose = []
            seps = []
            for dist,j in close:
                if dist < LIM:
                    jclose.append(j)
                    seps.append(dist)
                else:
                    break
            match.append(fromarrays([isorted[jclose], seps],
                                    dtype=[('ind','i8'),('sep','f8')]))

    if not allmatches:
        # return both indices and separations in a recarray
        temp = np.rec.fromrecords(match, names='ind,sep')
        # change to arcseconds
        temp.sep *= 3600.
        temp.sep[temp.sep < 0] = -1.
        return temp
    else:
        return match
