#!/bin/bash

file_list=$1
shift

cam=C0
filter=W
tag=
biasdir=
flatdir=

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ -s "$file_list" ] || { echo "No $file_list, exiting..." ; exit 1 ; }
[ "$biasdir" ] || { echo "No biasdir provided, exiting..." ; exit 1 ; }
[ "$flatdir" ] || { echo "No flatdir provided, exiting..." ; exit 1 ; }

go_iter=0
function gonogo() {
    ((go_iter++))
    [ "$((go_iter%NBATCH))" -eq 0 ] && wait
}

[ -d tmp ] || mkdir tmp

[ "$tag" ] && tag="${tag}_"
biasfile=${tag}bias_${cam}.fits
flatfile=${tag}flat_${cam}.fits

file0=`head -1 $file_list`
day=`basename $file0 | cut -c1-8`
here=`pwd`

if [ -f "$biasfile" ]; then
    echo "Using biasfile $biasfile"
else
    if [ -f ${biasdir}/$biasfile ]; then
        echo "Using biasfile ${biasdir}/$biasfile"
    else
        bias=`find_bias.sh $day $cam`
        cd $biasdir
        echo funpack -O bias_${cam}.fits $bias
        funpack -O bias_${cam}.fits $bias
        cd $here
    fi
    ln -s ${biasdir}/$biasfile .
fi

if [ -f "$flatfile" ]; then
    echo "Using flatfile $flatfile"
else
    if [ -f ${flatdir}/$flatfile ]; then
        echo "Using flatfile ${flatdir}/$flatfile"
    else
        flat=`find_flat.sh $day $cam $filter`
        cd $flatdir
        echo funpack -O flat_${cam}.fits $flat
        funpack -O flat_${cam}.fits $flat
        cd $here
    fi
    ln -s ${flatdir}/$flatfile .
fi

function doit() {
    local file=$1
    imreduce $file $biasfile $flatfile tmp/$file
    mv tmp/$file $file
}

ddoti_split.sh $file_list 1
for file in `cat $file_list`; do
    [ -f f00_$file ] && mv f00_$file $file
done

for file in `cat $file_list`; do
    doit $file &
    gonogo
done
wait
