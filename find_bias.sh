#!/bin/sh

day=$1
shift
ccd=$1
shift

dir=${REDUX_BASE_DIR}/bias_bank

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ -d "$dir" ] || mkdir $dir

[ "$ccd" = "C2" -o "$ccd" = "C3" ] && exit 1

cd $dir

ls bias*.fits* 2>/dev/null | grep $ccd > $dir/file_list$$.txt
cut -c6-13 file_list$$.txt > days$$.txt
paste days$$.txt file_list$$.txt > days$$.tmp
sort -n -k 1 days$$.tmp > days$$.txt
rm days$$.tmp

file=`awk '{if($1>'$day') exit; file=$2}END{print file}' days$$.txt`
echo ${dir}/$file
rm days$$.txt file_list$$.txt
