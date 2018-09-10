#!/bin/bash

file_list=$1
nsplit=$2

go_iter=0
function gonogo() {
    ((go_iter++))
    [ "$((go_iter%NBATCH))" -eq 0 ] && wait
}

ddoti_split.py $file_list $nsplit > split$$.batch

grep sethead split$$.batch > split$$.batch1
grep -v sethead split$$.batch > split$$.batch2

for n in `awk '{print NR}' split$$.batch1`; do
    cmd=`sed -n "${n}p" split$$.batch1`
    eval $cmd >/dev/null 2>&1 &
    gonogo
done
wait

for n in `awk '{print NR}' split$$.batch2`; do
    cmd=`sed -n "${n}p" split$$.batch2`
    eval $cmd >/dev/null 2>&1 &
    gonogo
done
wait

rm split$$.batch split$$.batch1 split$$.batch2
