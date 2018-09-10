#!/bin/bash

file_list=$1
shift

nstars_min=10
thresh=1.5
cal_mag_max=16.0
inst=nircam
sex_args=
srcfile=

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

go_iter=0
function gonogo() {
    ((go_iter++))
    [ "$((go_iter%NBATCH))" -eq 0 ] && wait
}

[ "$thresh" ] || thresh=1.5

if [ "$srcfile" ]; then
    for file in `cat $file_list`; do
        x=`gethead CRPIX1 $file`; y=`gethead CRPIX2 $file`
        base=${file%'.fits'}
        if [ -d ${base}_dir ]; then rm ${base}_dir/sky*
        else mkdir ${base}_dir ; fi
        awk '{print $1+'$x',$2+'$y',$3,$4,$5}' $srcfile > ${base}_dir/sky.list &
        gonogo
    done
    wait
fi

for file in `cat $file_list`; do
    run_sex.sh $file $inst -DETECT_THRESH $thresh $sex_args > /dev/null 2>&1 &
    gonogo
done
wait

if [ "$srcfile" ]; then

    for base in `sed -e 's/\.fits//g' $file_list`; do
        calibrate.py ${base}_dir/${base}_radec.txt $srcfile ${base}_dir/${base}_radec.txt.photometry.txt ${base}_dir/${base}_radec.txt.match.txt $cal_mag_max > ${base}_dir/${base}.report &
        gonogo
    done
    wait

    for file in `cat $file_list`; do
        base=${file%'.fits'}
        ns=`head -1 ${base}_dir/${base}.report | awk '{print $2}'`
        [ "$ns" ] || ns=0
        [ "$ns" -ge "$nstars_min" ] && echo $file
    done > stars_$file_list

else

    for file in `cat $file_list`; do
        base=${file%'.fits'}
        fwhm=`quick_mode ${base}_dir/${base}_radec.txt n=7`
        if [ "$fwhm" ]; then
            head -1 ${base}_dir/${base}_radec.txt > ${base}_dir/${base}_radec.tmp
            grep -v '#' ${base}_dir/${base}_radec.txt | awk '{if($7>'$fwhm'/2. && $7<'$fwhm'*4) print}' >> ${base}_dir/${base}_radec.tmp
            mv ${base}_dir/${base}_radec.tmp ${base}_dir/${base}_radec.txt
            nstars=`grep -v '#' ${base}_dir/${base}_radec.txt | wc -l`
            sethead NSTARS=$nstars FWHM=$fwhm $file
            echo $fwhm
        fi
    done > fwhm$$.txt
    fwhm0=`quick_mode fwhm$$.txt`
    [ "$fwhm0" ] && echo "Median FWHM of stars in images is $fwhm0 pixels"
    rm fwhm$$.txt

    gethead -a NSTARS @$file_list | awk '{if($2>'$nstars_min') print $1}' > stars_$file_list
fi

n0=`cat $file_list | wc -l`
n1=`cat stars_$file_list | wc -l`
echo "Keeping $n1 of $n0 images with $nstars_min or more stars."
