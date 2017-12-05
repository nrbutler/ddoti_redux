#!/bin/sh

day=$1
shift
ccd=$1
shift
filter=$1
shift

dir=${REDUX_BASE_DIR}/flat_bank

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ -d "$dir" ] || mkdir $dir

cd $dir

ls flat*.fits* 2>/dev/null > file_list0$$.txt
ls flat*.fits* 2>/dev/null | tr '_.' '  ' | awk '{print $3,$4}' > file_list1$$.txt
paste file_list0$$.txt file_list1$$.txt | eval "awk '{if(\$2~/$ccd/) print}'" > file_list$$.txt
if [ "$ccd" = "C0" ]; then
    eval "awk '{if(\$3~/$filter/) print}'" file_list$$.txt > file_list1$$.txt
    [ -s file_list1$$.txt ] || grep -i $filter file_list$$.txt > file_list1$$.txt
    mv file_list1$$.txt file_list$$.txt
fi

cut -c6-13 file_list$$.txt > days$$.txt
paste days$$.txt file_list$$.txt > days$$.tmp
sort -n -k 1 days$$.tmp > days$$.txt
rm days$$.tmp

file=`awk '{if($1>'$day') exit; file=$2}END{print file}' days$$.txt`
[ "$file" ] || file=`awk '{print $2; exit}' days$$.txt`
echo ${dir}/$file
rm days$$.txt file_list$$.txt file_list0$$.txt file_list1$$.txt 2>/dev/null
