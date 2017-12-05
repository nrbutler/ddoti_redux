#!/bin/bash

file_list=$1
shift

grb_RA=
grb_DEC=

#matching
match_rad=1000 # pixels, default 0 means use image size (slow)
match_thresh=20.0
mag_diff_max=3.0
nmatch_min=10
image_rad_max=1500

catdir=

# base reduction
sat_level=50000.0
gain=0.64
inst=nircam
do_distort=yes

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ -s $file_list ] || { echo "No file $file_list" ; exit 1 ; }
[ "$catdir" ] || { echo "No catalog dir specified" ; exit 2 ; }

go_iter=0
function gonogo() {
    ((go_iter++))
    [ "$((go_iter%NBATCH))" -eq 0 ] && wait
}

file0=`head -1 $file_list`
[ -f $file0 ] || { echo "Cannot find first file $file0" ; exit 3 ; }

file0f=`readlink -f $file0`
cam=`gethead CCD_NAME $file0 | cut -c1-2`
[ "$cam" ] || cam=`basename $file0 | cut -c16-17`
[ "$filter" ] || filter=`gethead FILTER $file0`
[ "$filter" ] || filter=W

sgain=`gethead SOFTGAIN $file0`
[ "$sgain" ] || sgain=1.0
gain=`echo $gain $sgain | awk '{print $1*$2}'`
sat_level=`echo $sat_level $sgain | awk '{print $1/$2}'`
gethead -a -u EXPTIME ETROBAM @$file_list | awk '{am=$3;if(am~/___/) am=1.0; printf("sethead EXPTIME0=%.2f AIRMASS=%.2f AM_COEF=0.1 %s\n",$2,am,$1)}' | sh
sethead SATURATE=$sat_level GAIN=$gain @$file_list

gethead -a -u CTYPE1 NAXIS1 @$file_list | eval awk "'{if(\$2~/___/) print \$1>\"new_$file_list\"; else print \$1>\"old_$file_list\"}'"
nnew=`cat new_$file_list 2>/dev/null | wc -l`
[ "$nnew" -eq 0 ] && { echo "No new files to align, exiting..." ; exit 4 ; }

#
# starting work
#

sethead CTYPE1='RA---TAN' CTYPE2='DEC--TAN' @new_$file_list

base0=${file0%'.fits'}
if [ -f ${base0}.wcs ]; then
    echo "Already ran astrometry.net on $file0"
    do_distort=no
    cpkeys="CRVAL1 CRVAL2 CD1_1 CD1_2 CD2_1 CD2_2 CRPIX1 CRPIX2"
    cpkeys="$cpkeys PV1_1 PV2_1 PV1_7 PV1_9 PV2_7 PV2_9 PV1_17 PV1_21 PV2_17 PV2_21 PV1_19 PV2_19"
    cphead $cpkeys $file0 @new_$file_list
else
    run_astnet_blind.sh $file0 ra=$grb_RA dec=$grb_DEC inst=$inst detect_thresh=$match_thresh
    [ -f ${base0}.wcs ] || { echo "Astrometry.net failed for $file0" ; exit 5 ; }

    cphead CRVAL1 CRVAL2 CRPIX1 CRPIX2 CD1_1 CD1_2 CD2_1 CD2_2 ${base0}.wcs @new_$file_list
    sethead PV1_1=1.0 PV2_1=1.0 PV1_7=3.3E-4 PV1_9=3.3E-4 PV2_7=3.3E-4 PV2_9=3.3E-4 PV1_17=2.0E-5 PV1_21=2.0E-5 PV2_17=2.0E-5 PV2_21=2.0E-5 PV1_19=4.0E-5 PV2_19=4.0E-5 @new_$file_list
fi
grb_RA=`gethead CRVAL1 $file0 | awk '{printf("%.8f\n",$1)}'`
grb_DEC=`gethead CRVAL2 $file0 | awk '{printf("%.8f\n",$1)}'`

find_stars.sh new_$file_list inst=$inst thresh=$match_thresh nstars_min=10 &
pid0=$!

if [ -f "${catdir}/usno_radec16_${cam}.txt" ]; then
    echo "Using ${catdir}/usno_radec16_${cam}.txt"
else
    here=`pwd`
    cd $catdir
    echo grab_usno_local.sh $file0f
    grab_usno_local.sh $file0f
    cd $here
fi
ln -s ${catdir}/usno_radec16_${cam}.txt usno_radec16.txt 2>/dev/null

wait $pid0
mv stars_new_$file_list new_$file_list

match_ref=usno_radec16.txt
for file in `cat new_$file_list`; do
    base=${file%'.fits'}
    sexlist_cc.py $file ${base}_dir/${base}_radec.txt $match_ref X $match_rad 0 0 $mag_diff_max $nmatch_min $image_rad_max >> match$$.report &
    gonogo
done
wait
cat match$$.report

awk '{print $6}' match$$.report | awk -F/ '{print $2}' | sed -e 's/_radec\.txt/\.fits/g' | sort > new_$file_list
cat old_$file_list new_$file_list 2>/dev/null > $file_list

[ "$do_distort" = "yes" ] && ddoti_distortion.py $file_list

rm match$$.report
