#!/bin/sh

file=$1
shift

overfac=1.5

cat_dir=${REDUX_BASE_DIR}/catalogs/usno

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

cam=`gethead CCD_NAME $file | cut -c1-2`
[ "$cam" ] || cam=`basename $file | cut -c16-17`
[ "$cam" ] || cam=X

ra=`gethead CRVAL1 $file`
dec=`gethead CRVAL2 $file`

dx=`gethead NAXIS1 $file`
dy=`gethead NAXIS2 $file`

cdec=`echo $dec | awk '{print cos($1/57.296)}'`
pscale=`gethead CD1_1 CD1_2 $file | awk '{print sqrt($1*$1+$2*$2)}'`
dra=`echo $dx $pscale $overfac $cdec | awk '{print $1*$2*$3/$4}'`
ddec=`echo $dy $pscale $overfac | awk '{print $1*$2*$3}'`

ras=`echo "$ra $dra"  | awk '{r0=$1-$2/2+360; r1=$1+$2/2+360; r00 = int(r0/10)*10+5; r11 = int(r1/10)*10+5; for (i=r00;i<=r11;i+=10) if (i<360) print i; else print i-360}'`
decs=`echo "$dec $ddec"  | awk '{d0=$1-$2/2+180; d1=$1+$2/2+180; d00 = int(d0/2)*2+1; d11 = int(d1/2)*2+1; for(i=d00;i<=d11;i+=2) if(i<90) print i; else print i-180}'`

for ra0 in $ras ; do
    for dec0 in $decs ; do
        echo grab_usno_local.py ${cat_dir}/usno_b1_${ra0}_${dec0}.fits $ra $dra $dec $ddec
        grab_usno_local.py ${cat_dir}/usno_b1_${ra0}_${dec0}.fits $ra $dra $dec $ddec > udata$$_${ra0}_${dec0}.txt &
    done
done
wait

cat udata$$_*.txt > usno_radec_${cam}.txt
rm udata$$_*.txt

wc -l usno_radec_${cam}.txt
awk '{if($3<=16) print}' usno_radec_${cam}.txt > usno_radec16_${cam}.txt
