#!/bin/bash

file_list=$1
shift

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

file0=`head -1 $file_list`
cam=`basename $file0 | cut -c16-17`
bin=`gethead SDTBN $file0`
[ "$bin" ] || bin=1
day=`basename $file0 | cut -c1-8`
filter=`gethead FILTER $file0`
[ "$filter" ] || filter=W
here=`pwd`

biasfile=bias_${cam}.fits
flatfile=flat_${cam}.fits

if [ -f "$biasfile" ]; then
    echo "Using biasfile $biasfile"
else
    if [ -f ${biasdir}/$biasfile ]; then
        echo "Using biasfile ${biasdir}/$biasfile"
    else
        bias=`find_bias.sh $day $cam $bin`
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
        flat=`find_flat.sh $day $cam $filter $bin`
        cd $flatdir
        echo funpack -O flat_${cam}.fits $flat
        funpack -O flat_${cam}.fits $flat
        cd $here
    fi
    ln -s ${flatdir}/$flatfile .
fi

# first pull out data_sec, and estimate bias from bias_sec (creates f00_file for each file)
ddoti_split.sh $file_list 1

for file in `cat $file_list`; do
    imreduce f00_$file $biasfile $flatfile \!$file &
    gonogo
done
wait

rm f00_*.fits 2>/dev/null
