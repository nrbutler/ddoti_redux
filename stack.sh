#!/bin/sh

file_list=$1
shift
outfile=$1
shift

combine_type=MEDIAN

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

echo swarp @$file_list -c ${SWARP_DIR}/ddoti_redux.swarp -COMBINE_TYPE $combine_type -IMAGEOUT_NAME $outfile -WEIGHTOUT_NAME expmap$outfile -WEIGHT_TYPE NONE -RESAMPLE N -SUBTRACT_BACK N
swarp @$file_list -c ${SWARP_DIR}/ddoti_redux.swarp -COMBINE_TYPE $combine_type -IMAGEOUT_NAME $outfile -WEIGHTOUT_NAME expmap$outfile -WEIGHT_TYPE NONE -RESAMPLE N -SUBTRACT_BACK N 2>/dev/null
