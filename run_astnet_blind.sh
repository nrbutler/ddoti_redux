#!/bin/sh

file=$1
shift

radius=  # in degrees, leave blank for blind search
inst=nircam
detect_thresh=100.0
nstars=1000
ra=
dec=

cleanup=yes

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ -f $file ] || { echo "No file $file" ; exit 1 ; }

backend=${ASTNET_DIR}/ddoti_astnet_backend.cfg

[ "$ra" ] || ra=`gethead CRVAL1 $file`
[ "$dec" ] || dec=`gethead CRVAL2 $file`

search=
[ "$radius" ] && search="--ra $ra --dec $dec --radius $radius"

pscale=`gethead CD1_1 CD1_2 $file | awk '{print sqrt($1*$1+$2*$2)*3600.}'`
pscale1=`echo $pscale | awk '{print 0.9*$1}'`
pscale2=`echo $pscale | awk '{print 1.1*$1}'`
scale="--scale-units arcsecperpix --scale-low $pscale1 --scale-high $pscale2"

nx=`gethead NAXIS1 $file`
ny=`gethead NAXIS2 $file`
x0=`gethead CRPIX1 $file`
[ "$x0" ] || x0=`echo $nx | awk '{printf("%.0f\n",$1/2.)}'`
y0=`gethead CRPIX2 $file`
[ "$y0" ] || y0=`echo $ny | awk '{printf("%.0f\n",$1/2.)}'`
crpix="--crpix-x $x0 --crpix-y $y0"

base=${file%'.fits'}
[ -d ${base}_dir ] || run_sex.sh $file $inst -DETECT_THRESH $detect_thresh
grep -v '#' ${base}_dir/${base}_radec.txt | sort -n -k 3 | awk '{if(NR<='$nstars') print $8,$9}' > xy$$.txt

[ -f ${base}.wcs ] && rm ${base}.wcs

echo 'X E "" "" "" "" "" ""
Y E "" "" "" "" "" ""' > hf$$.txt

echo "import pyfits
x=pyfits.tableload('xy$$.txt','hf$$.txt')
x.writeto('${base}.xy',clobber=True)" | python

echo "solve-field --backend-config $backend ${base}.xy --continue -w $nx -e $ny $crpix $scale $search --no-verify --no-plots --cpulimit 30 -T --pnm out.pnm --new-fits none"
solve-field --backend-config $backend ${base}.xy --continue -w $nx -e $ny $crpix $scale $search --no-verify --no-plots --cpulimit 30 -T --pnm out.pnm --new-fits none 2>&1

if [ "$cleanup" = "yes" ]; then
    rm ${base}.xy ${base}.axy xy$$.txt hf$$.txt ${base}.rdls ${base}.match ${base}.corr ${base}-indx.xyls ${base}.solved
fi
