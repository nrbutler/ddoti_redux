#!/bin/bash

file_list=$1
shift

grb_RA=
grb_DEC=

#matching
match_thresh=20.0
max_offset=2.0 # degrees

catdir=

# base reduction
inst=nircam
do_distort=yes

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ -s $file_list ] || { echo "No file $file_list" ; exit 1 ; }
[ "$catdir" ] || { echo "No catalog dir specified" ; exit 2 ; }

file0=`head -1 $file_list`
[ -f $file0 ] || { echo "Cannot find first file $file0" ; exit 3 ; }

go_iter=0
function gonogo() {
    ((go_iter++))
    [ "$((go_iter%NBATCH))" -eq 0 ] && wait
}

cam=`gethead CCD_NAME $file0 | cut -c1-2`
[ "$cam" ] || cam=`basename $file0 | cut -c16-17`
[ "$filter" ] || filter=`gethead FILTER $file0`
[ "$filter" ] || filter=W

gethead -a -u EXPTIME ETROBAM @$file_list | awk '{am=$3;if(am~/___/) am=1.0; printf("sethead EXPTIME0=%.2f AIRMASS=%.2f AM_COEF=0.1 %s\n",$2,am,$1)}' | sh

gethead -a -u CTYPE1 NAXIS1 @$file_list | eval awk "'{if(\$2~/___/) print \$1>\"new_$file_list\"; else print \$1>\"old_$file_list\"}'"
# ignore images with pier flips during the exposure
gethead -a SMTMRO EMTMRO @new_$file_list > rot$$.txt
awk '{if($2==$3) print $1}' rot$$.txt > new_$file_list
rm rot$$.txt

nnew=`cat new_$file_list 2>/dev/null | wc -l`
[ "$nnew" -eq 0 ] && { echo "No new files to align, exiting..." ; exit 4 ; }

#
# starting work
#

echo find_stars.sh new_$file_list inst=$inst thresh=$match_thresh nstars_min=100
find_stars.sh new_$file_list inst=$inst thresh=$match_thresh nstars_min=100
mv stars_new_$file_list new_$file_list

# find the best file
gethead -a FWHM NSTARS @new_$file_list | awk '{print $1,$2*(3000./$3)}' | sort -n -k 2 | awk '{print $1}' > tmp_$file_list
mv tmp_$file_list new_$file_list
file0=`head -1 new_$file_list`

p1=4.6E-4   # expected barrel distortion 13*(rad/2)^3 arcsec (corrects 13 arcsec shrinkage at 2 degrees field angle)
ra0=$grb_RA ; dec0=$grb_DEC
if [ -f best.wcs ]; then
   echo "Using base alignment from best.wcs"
else
    run_astnet_blind.sh $file0 ra=$grb_RA dec=$grb_DEC radius=$max_offset p1=$p1 catdir=$catdir do_distort=$do_distort
    wfile0=${file0%'.fits'}.wcs
    [ -f "$wfile0" ] || { echo "Astrometry.net failed for $file0" ; rm new_$file_list 2>/dev/null ; exit 5 ; }
    ln -s $wfile0 best.wcs
fi

p1=`gethead PV1_7 best.wcs`
ra0=`gethead CRVAL1 best.wcs`; dec0=`gethead CRVAL2 best.wcs`
x0=`gethead CRPIX1 best.wcs`; y0=`gethead CRPIX2 best.wcs`
parity=`gethead CD1_1 CD2_2 CD1_2 CD2_2 best.wcs | awk '{p=$1*$2-$3*$4; if (p<0) print "pos"; else print "neg"}'`
for file in `cat new_$file_list`; do
    base=${file%'.fits'}
    [ -f "${base}.wcs" ] || {
        run_astnet_blind.sh $file ra=$ra0 dec=$dec0 x0=$x0 y0=$y0 radius=$max_offset parity=$parity p1=$p1 catdir=$catdir &
        gonogo
    }
done
wait

for file in `cat new_$file_list`; do
    [ -f "${file%'.fits'}.wcs" ] && echo $file
done > tmp_new_$file_list
mv tmp_new_$file_list new_$file_list

cat old_$file_list new_$file_list 2>/dev/null > $file_list
