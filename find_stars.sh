#!/bin/bash

file_list=$1
shift

nstars_min=1
thresh=1.5
inst=nircam
psf_only=no
plot=no
dao=
srcfile=

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

go_iter=0
function gonogo() {
    ((go_iter++))
    [ "$((go_iter%NBATCH))" -eq 0 ] && wait
}

[ "$thresh" ] || thresh=1.5

if [ "$srcfile" ]; then
    srcfile0=$srcfile
    srcfile="srcfile=$srcfile"
fi

if [ "$plot" = "yes" ]; then
    plot="_plot"
else
    plot=
fi

psfonly=
if [ "$psf_only" = "yes" ]; then
    psfonly=psfonly
fi

n0=`cat $file_list | wc -l`

function dophot() {
    local file=$1
    local base=${file%'.fits'}
    [ -d ${base}_dir ] && rm -r ${base}_dir
    run_sex.sh $file $inst $dao $psfonly $srcfile -DETECT_THRESH $thresh > /dev/null 2>&1
    [ "$srcfile" ] && calibrate.py ${base}_dir/${base}_radec.txt $srcfile0 > ${base}_dir/${base}.report
}

for file in `cat $file_list`; do
    dophot $file &
    gonogo
done
wait

if [ "$srcfile" ]; then

    for file in `cat $file_list`; do
        base=${file%'.fits'}
        ns=`head -1 ${base}_dir/${base}.report | awk '{print $2}'`
        [ "$ns" ] || ns=0
        [ "$ns" -ge "$nstars_min" ] && echo $file
    done > stars_$file_list

else

    for file in `cat $file_list`; do
        base=${file%'.fits'}
        grep -v '#' ${base}_dir/${base}_radec.txt | awk '{print $7}'
    done > fwhm$$.txt
    fwhm0=`quick_mode fwhm$$.txt`
    [ "$fwhm0" ] && echo "Median FWHM of stars in images is $fwhm0 pixels"
    rm fwhm$$.txt

    for file in `cat $file_list`; do
        base=${file%'.fits'}
        head -1 ${base}_dir/${base}_radec.txt > ${base}_dir/${base}_radec.tmp
        grep -v '#' ${base}_dir/${base}_radec.txt | awk '{if($7>'$fwhm0'/2. && $7<'$fwhm0'*4) print}' >> ${base}_dir/${base}_radec.tmp
        mv ${base}_dir/${base}_radec.tmp ${base}_dir/${base}_radec.txt
        ns=`grep -v '#' ${base}_dir/${base}_radec.txt | wc -l`
        [ "$ns" ] || ns=0
        [ "$ns" -ge "$nstars_min" ] && echo $file
    done > stars_$file_list

fi

n1=`cat stars_$file_list | wc -l`
echo "Keeping $n1 of $n0 images with $nstars_min or more stars."
