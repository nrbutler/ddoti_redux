#!/bin/bash

file_list=$1
shift

# base reduction
nsplit=4
wt_file=

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ -s $file_list ] || { echo "No file $file_list" ; exit 1 ; }
[ "$wt_file" ] || { echo "Must specify a weight file wt_file" ; exit 2 ; }
[ -f "$wt_file" ] || { echo "Weight file $wt_file does not exist" ; exit 3 ; }

go_iter=0
function gonogo() {
    ((go_iter++))
    [ "$((go_iter%NBATCH))" -eq 0 ] && wait
}

file0=`head -1 $file_list | sed -e 's/\.bz2/\.txt/g'`
[ -f $file0 ] || { echo "Cannot find first file $file0" ; exit 3 ; }
cam=`gethead CCD_NAME $file0 | cut -c1-2`
[ "$cam" ] || cam=`basename $file0 | cut -c16-17`

sargs="-c ${SWARP_DIR}/ddoti_redux.swarp -RESAMPLING_TYPE NEAREST -BACK_SIZE 256 -COMBINE N"
sargs="$sargs -IMAGEOUT_NAME stack_${cam}.fits -WEIGHTOUT_NAME stack_${cam}.wt.fits"
rm stack_${cam}.head 2>/dev/null
echo swarp @$file_list $sargs -HEADER_ONLY Y -WEIGHT_IMAGE $wt_file
swarp @$file_list $sargs -HEADER_ONLY Y -WEIGHT_IMAGE $wt_file

if [ -f stack_${cam}.head0 ]; then
    echo "Using existing projection from stack_${cam}.head0"

    ps=`gethead CD1_1 CD1_2 stack_${cam}.head0 | awk '{print sqrt($1*$1+$2*$2)}'`
    ra0=`gethead CRVAL1 stack_${cam}.head0`; dec0=`gethead CRVAL2 stack_${cam}.head0`
    x0=`gethead CRPIX1 stack_${cam}.head0`; y0=`gethead CRPIX2 stack_${cam}.head0`
    ra1=`gethead CRVAL1 stack_${cam}.fits`; dec1=`gethead CRVAL2 stack_${cam}.fits`
    x1=`gethead CRPIX1 stack_${cam}.fits`; y1=`gethead CRPIX2 stack_${cam}.fits`
    cdec0=`echo $dec0 | awk '{print cos($1/57.3)}'`

    dx=`echo $ra1 $ra0 $x1 $x0 | awk '{dx=($1-$2)*'$cdec0'/'$ps'-($4-$3); if(dx>0) printf("%.0f\n",dx); else print 0}'`
    dy=`echo $dec1 $dec0 $y1 $y0 | awk '{dx=($2-$1)/'$ps'-($4-$3); if(dx>0) printf("%.0f\n",dx); else print 0}'`

    nx=`gethead NAXIS1 stack_${cam}.fits`
    ny=`gethead NAXIS2 stack_${cam}.fits`
    x0=`gethead CRPIX1 stack_${cam}.head0 | awk '{print $1+'$dx'}'`
    y0=`gethead CRPIX2 stack_${cam}.head0 | awk '{print $1+'$dy'}'`
    cp stack_${cam}.head0 stack_${cam}.head
    sethead NAXIS1=$nx NAXIS2=$ny CRPIX1=$x0 CRPIX2=$y0 stack_${cam}.head
else
    imhead -f stack_${cam}.fits > stack_${cam}.head
    cp stack_${cam}.head stack_${cam}.head0
fi

for file in `cat new_$file_list`; do
    echo swarp $file $sargs -WEIGHT_IMAGE $wt_file
    swarp $file $sargs -WEIGHT_IMAGE $wt_file 2>/dev/null &
    gonogo
done
wait

ddoti_split.py stack_${cam}.head $nsplit 2>/dev/null | sh

sargs="-c ${SWARP_DIR}/ddoti_redux.swarp -COMBINE_TYPE MEDIAN -WEIGHT_SUFFIX .weight.fits -RESAMPLE N"
sed -e 's/\.fits/\.resamp\.fits/' $file_list > r$file_list
for tag0 in `ls f??_stack_${cam}.head | awk -F_ '{print $1}'`; do
    echo swarp @r$file_list $sargs -IMAGEOUT_NAME ${tag0}_stack_${cam}.fits -WEIGHTOUT_NAME ${tag0}_stack_${cam}.wt.fits
    swarp @r$file_list $sargs -IMAGEOUT_NAME ${tag0}_stack_${cam}.fits -WEIGHTOUT_NAME ${tag0}_stack_${cam}.wt.fits 2>/dev/null &
    gonogo
done
wait

# finally, combine the stacks into a master stack
ls f??_stack_${cam}.fits > stack_$file_list
sargs="-c ${SWARP_DIR}/ddoti_redux.swarp -COMBINE_TYPE AVERAGE -WEIGHT_SUFFIX .wt.fits -RESAMPLE N"
echo swarp @stack_$file_list $sargs -IMAGEOUT_NAME stack_${cam}.fits -WEIGHTOUT_NAME stack_${cam}.wt.fits
swarp @stack_$file_list $sargs -IMAGEOUT_NAME stack_${cam}.fits -WEIGHTOUT_NAME stack_${cam}.wt.fits 2>/dev/null

exptime=`gethead -a EXPTIME @$file_list | awk '{s=s+$2}END{printf("%.2f\n",s)}'`
sat_level=`gethead SATURATE $file0`
sethead EXPTIME=$exptime SATURATE=$sat_level stack_${cam}.fits stack_${cam}.wt.fits
sethead EXPTIME=$exptime SATURATE=$sat_level @stack_$file_list

rm f??_stack_${cam}.head rf??_$file_list r$file_list 2>/dev/null
