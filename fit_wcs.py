from numpy import sqrt,sin,cos,sqrt,arctan2,arcsin,arccos

def xy2ij(x,y,hdr,distort=False):

    cd1_1, cd1_2 = hdr['CD1_1'], hdr['CD1_2']
    cd2_1, cd2_2 = hdr['CD2_1'], hdr['CD2_2']
    x0, y0 = hdr['CRPIX1'], hdr['CRPIX2']

    i = cd1_1*(x-x0)+cd1_2*(y-y0)
    j = cd2_1*(x-x0)+cd2_2*(y-y0)

    if (distort):

        try: a = hdr['PV1_7']
        except: a=0.

        fac=a*(i**2+j**2)
        i += fac*i
        j += fac*j

    return i,j


def ij2ad(i,j,hdr):

    fac = 57.29577951

    xi,eta = i/fac,j/fac

    xi0, eta0 = hdr['CRVAL1']/fac, hdr['CRVAL2']/fac
    seta0, ceta0 = sin(eta0), cos(eta0)

    cosc = 1./sqrt(1.+xi**2.+eta**2.)
    a = fac*( xi0 + arctan2( xi, ceta0 - eta*seta0 ) )
    d = fac*( arcsin(cosc*(seta0 + eta*ceta0)) )

    try:
        if (a<0): a+=360
        elif (a>360): a-=360
    except:
        a[a<0]+=360
        a[a>360]-=360

    return a,d


def ad2ij(a,d,hdr):
    fac = 57.29577951

    a0, d0 = hdr['CRVAL1']/fac, hdr['CRVAL2']/fac
    sa,ca = sin(a/fac),cos(a/fac)
    sd,cd = sin(d/fac),cos(d/fac)
    sa0,ca0 = sin(a0),cos(a0)
    sd0,cd0 = sin(d0),cos(d0)

    cos_c = sd0*sd + cd0*cd*(ca*ca0+sa*sa0)
    i = fac*cd*(sa*ca0-ca*sa0) / cos_c
    j = fac*(cd0*sd-sd0*cd*(ca*ca0+sa*sa0)) / cos_c

    return i,j


def ij2xy(i,j,hdr,distort=False):

    cd1_1, cd1_2 = hdr['CD1_1'], hdr['CD1_2']
    cd2_1, cd2_2 = hdr['CD2_1'], hdr['CD2_2']
    denom = cd1_1*cd2_2 - cd1_2*cd2_1

    x0, y0 = hdr['CRPIX1'], hdr['CRPIX2']

    if (distort):

        try: a = hdr['PV1_7']
        except: a=0.

        fac=a*(i**2+j**2)
        i1,j1 = i/(1.+fac),j/(1.+fac)
        for ii in range(3):
            fac=a*(i1**2+j1**2)
            i1,j1 = i/(1.+fac),j/(1.+fac)

        x = x0 + (i1*cd2_2 - j1*cd1_2)/denom
        y = y0 + (j1*cd1_1 - i1*cd2_1)/denom

    else:

        x = x0 + (i*cd2_2 - j*cd1_2)/denom
        y = y0 + (j*cd1_1 - i*cd2_1)/denom

    return x,y


def xy2ad(x,y,hdr,distort=False):
    i,j = xy2ij(x,y,hdr,distort=distort)
    return ij2ad(i,j,hdr)


def ad2xy(a,d,hdr,distort=False):
    i,j = ad2ij(a,d,hdr)
    return ij2xy(i,j,hdr,distort=distort)
