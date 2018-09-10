#!/bin/sh

file=$1
shift

radius=2.0  # in degrees
inst=nircam
detect_thresh=20.0
nstars=1000
# leave these blank for blind search:
ra=
dec=
x0=
y0=
parity=

# distortion
p1=
do_distort=no

# for wcs refinement
catdir=.

cleanup=yes

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ -f $file ] || { echo "No file $file" ; exit 1 ; }

cam=`gethead CCD_NAME $file | cut -c1-2`
[ "$cam" ] || cam=`basename $file | cut -c16-17`

backend=${ASTNET_DIR}/ddoti_astnet_backend.cfg

[ "$parity" ] && parity="--parity $parity"

search=
[ "$ra" -a "$dec" ] && search="--ra $ra --dec $dec --radius $radius"

bin=`gethead SDTBN $file`
pscale1=`echo $bin | awk '{print $1*1.9}'`
pscale2=`echo $bin | awk '{print $1*2.1}'`
scale="--scale-units arcsecperpix --scale-low $pscale1 --scale-high $pscale2"

nx=`gethead NAXIS1 $file`
ny=`gethead NAXIS2 $file`
[ "$x0" ] || x0=`echo $nx | awk '{printf("%.0f\n",$1/2.)}'`
[ "$y0" ] || y0=`echo $ny | awk '{printf("%.0f\n",$1/2.)}'`
crpix="--crpix-x $x0 --crpix-y $y0"

base=${file%'.fits'}
[ -d ${base}_dir ] || run_sex.sh $file $inst -DETECT_THRESH $detect_thresh
grep -v '#' ${base}_dir/${base}_radec.txt | sort -n -k 3 | awk '{if(NR<='$nstars') print $8,$9}' > xy$$.txt

[ -f ${base}.wcs ] && rm ${base}.wcs

echo 'X E "" "" "" "" "" ""
Y E "" "" "" "" "" ""' > hf$$.txt

echo "from pyfits import tableload
x=tableload('xy$$.txt','hf$$.txt')
x.writeto('${base}.xy',clobber=True)" | python

echo "solve-field --backend-config $backend ${base}.xy --continue -w $nx -e $ny $crpix $scale $search --no-verify --no-plots --cpulimit 3 -T --pnm out.pnm --new-fits none $parity"
solve-field --backend-config $backend ${base}.xy --continue -w $nx -e $ny $crpix $scale $search --no-verify --no-plots --cpulimit 3 -T --pnm out.pnm --new-fits none $parity 2>&1

if [ -f "${base}.wcs" ]; then
    # now go in and update ra and dec values in ${base}_dir/${base}_radec.txt using ${base}.wcs
    sethead CTYPE1='RA---TAN' CTYPE2='DEC--TAN' $file
    cphead CRVAL1 CRVAL2 CD1_1 CD1_2 CD2_1 CD2_2 CRPIX1 CRPIX2 ${base}.wcs $file
    [ "$p1" ] && sethead PV1_1=1.0 PV2_1=1.0 PV1_7=$p1 PV1_9=$p1 PV2_7=$p1 PV2_9=$p1 $file

    if [ -f "${catdir}/usno_radec16_${cam}.fits" ]; then
        echo "Using ${catdir}/usno_radec16_${cam}.fits"
    else
        here=`pwd`
        filef=`readlink -f $file`
        cd $catdir
        echo grab_usno_local.sh $filef max_offset=$radius
        grab_usno_local.sh $filef max_offset=$radius
        cd $here
    fi
    ln -s ${catdir}/usno_radec_${cam}.fits usno_radec.fits 2>/dev/null
    ln -s ${catdir}/usno_radec_${cam}.pm.fits usno_radec.pm.fits 2>/dev/null
    ln -s ${catdir}/usno_radec16_${cam}.fits usno_radec16.fits 2>/dev/null

    echo tune_wcs.py $file usno_radec16.fits
    tune_wcs.py $file usno_radec16.fits

    if [ "$do_distort" = "yes" ]; then
        echo ddoti_distortion.py $file $p1
        ddoti_distortion.py $file $p1
        p1=`gethead PV1_7 $file`
    fi
    x0=`gethead CRPIX1 $file`; y0=`gethead CRPIX2 $file`
    ra0=`gethead CRVAL1 $file`; dec0=`gethead CRVAL2 $file`
    sethead PV1_7=$p1 CRPIX1=$x0 CRPIX2=$y0 CRVAL1=$ra0 CRVAL2=$dec0 ${base}.wcs
fi

if [ "$cleanup" = "yes" ]; then
    rm ${base}.xy ${base}.axy xy$$.txt hf$$.txt ${base}.rdls ${base}.match ${base}.corr ${base}-indx.xyls ${base}.solved 2>/dev/null
fi
