#!/bin/bash

file_list=$1
shift

# base reduction
nsplit=4
wt_file=
inst=nircam
mask_sigma=3.0

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ -s $file_list ] || { echo "No file $file_list" ; exit 1 ; }
[ "$wt_file" ] || { echo "Must specify a weight file wt_file" ; exit 2 ; }
[ -f "$wt_file" ] || { echo "Weight file $wt_file does not exist" ; exit 3 ; }

file0=`head -1 $file_list | sed -e 's/\.bz2/\.txt/g' -e 's/\.fz/\.txt/g'`
[ -f $file0 ] || { echo "Cannot find first file $file0" ; exit 3 ; }

go_iter=0
function gonogo() {
    ((go_iter++))
    [ "$((go_iter%NBATCH))" -eq 0 ] && wait
}

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

# use only the best files
n0=`cat $file_list | wc -l`
gethead -a FWHM NSTARS @$file_list | awk '{print $1,$2*(3000./$3)}' | sort -n -k 2 > file_stats$$.txt
stat0=`head -1 file_stats$$.txt | awk '{print $2}'`
awk '{if($2<2*'$stat0') print $1}' file_stats$$.txt > good_$file_list
n1=`cat good_$file_list | wc -l`
echo "Stacking $n1 of $n0 files with FWHM*(3000/Nstars) < 2*$stat0"

sargs="-c ${SWARP_DIR}/ddoti_redux.swarp -COMBINE_TYPE MEDIAN -WEIGHT_SUFFIX .weight.fits -RESAMPLE N"
ls `sed -e 's/\.fits/\.resamp\.fits/' good_$file_list` 2>/dev/null > r$file_list
for tag0 in `ls f??_stack_${cam}.head | awk -F_ '{print $1}'`; do
    echo swarp @r$file_list $sargs -IMAGEOUT_NAME ${tag0}_stack_${cam}.fits -WEIGHTOUT_NAME ${tag0}_stack_${cam}.wt.fits
    swarp @r$file_list $sargs -IMAGEOUT_NAME ${tag0}_stack_${cam}.fits -WEIGHTOUT_NAME ${tag0}_stack_${cam}.wt.fits 2>/dev/null &
    gonogo
done
wait

# combine the stacks into a master stack
ddoti_unsplit.py stack_${cam}.fits

exptime=`gethead -a EXPTIME @$file_list | awk '{s=s+$2}END{printf("%.2f\n",s)}'`
skylev=`gethead -a SKYLEV @$file_list | awk '{n=n+1;s=s+1/$2}END{printf("%.2f\n",n/s)}'`
sat_level=`gethead -a SATURATE @$file_list | awk '{n=n+1;s=s+$2}END{printf("%.2f\n",s/n)}'`
sethead EXPTIME=$exptime SKYLEV=$skylev SATURATE=$sat_level stack_${cam}.fits stack_${cam}.wt.fits
sethead EXPTIME=$exptime SKYLEV=$skylev SATURATE=$sat_level f??_stack_${cam}.fits f??_stack_${cam}.wt.fits

# make the source mask
[ -d stack_${cam}_dir ] && rm -r stack_${cam}_dir
run_sex.sh stack_${cam}.fits $inst -DETECT_THRESH 10.0 -CHECKIMAGE_TYPE OBJECTS -CHECKIMAGE_NAME mask.fits
sources2mask.py stack_${cam}.wt.fits stack_${cam}_dir/mask.fits maskstack_${cam}.fits $mask_sigma

x0=`gethead CRPIX1 stack_${cam}.fits`; nx0=`gethead NAXIS1 stack_${cam}.fits`
y0=`gethead CRPIX2 stack_${cam}.fits`; ny0=`gethead NAXIS2 stack_${cam}.fits`

function do_weights() {
    local file=$1
    local base=${file%'.fits'}
    local x=`gethead CRPIX1 $file`
    local y=`gethead CRPIX2 $file`
    local dx=`echo $x | awk '{printf("%.0f\n",'$x0'-$1+1)}'`
    local dy=`echo $y | awk '{printf("%.0f\n",'$y0'-$1+1)}'`
    local x1=`gethead NAXIS1 $file | awk '{printf("%.0f\n",'$dx'+$1-1)}'`
    local y1=`gethead NAXIS2 $file | awk '{printf("%.0f\n",'$dy'+$1-1)}'`
    [ "$x1" -gt "$nx0" ] && {
        dx=`echo $dx | awk '{printf("%.0f\n",$1-'$x1'+'$nx0')}'`
        x1=$nx0
    }
    [ "$y1" -gt "$ny0" ] && {
        dy=`echo $dy | awk '{printf("%.0f\n",$1-'$y1'+'$ny0')}'`
        y1=$ny0
    }
    echo weight_clip $file ${base}.weight.fits maskstack_${cam}.fits[$dx:$x1,$dy:$y1] stack_${cam}.fits[$dx:$x1,$dy:$y1] stack_${cam}.wt.fits[$dx:$x1,$dy:$y1]
    weight_clip $file ${base}.weight.fits maskstack_${cam}.fits[$dx:$x1,$dy:$y1] stack_${cam}.fits[$dx:$x1,$dy:$y1] stack_${cam}.wt.fits[$dx:$x1,$dy:$y1]
}

for file in `cat r${file_list}`; do
    do_weights $file &
    gonogo
done
wait

cp stack_${cam}.fits stack1_${cam}.fits
cp stack_${cam}.wt.fits stack1_${cam}.wt.fits

for tag0 in `ls f??_stack_${cam}.head | awk -F_ '{print $1}'`; do
    echo swarp @r$file_list $sargs -IMAGEOUT_NAME ${tag0}_stack_${cam}.fits -WEIGHTOUT_NAME ${tag0}_stack_${cam}.wt.fits
    swarp @r$file_list $sargs -IMAGEOUT_NAME ${tag0}_stack_${cam}.fits -WEIGHTOUT_NAME ${tag0}_stack_${cam}.wt.fits 2>/dev/null &
    gonogo
done
wait

# combine the stacks into a master stack
ddoti_unsplit.py stack_${cam}.fits

rm f??_stack_${cam}.head rf??_$file_list r$file_list 2>/dev/null
