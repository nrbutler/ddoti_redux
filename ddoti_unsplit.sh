#!/bin/bash

file_list=$1
nsplit=$2

n=0
for file in `cat $file_list`; do
    ddoti_unsplit.py $file $nsplit &
    ((n++))
    [ "$((n%NBATCH))" -eq 0 ] && wait
done
wait
