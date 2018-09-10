#!/bin/bash

file_list=$1
shift

catdir=
cam=C0

# source detection
do_lightcurves=yes
ndetect_max=50
search_snr=10
match_file=usno_radec.fits

# base reduction
inst=nircam

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ -s $file_list ] || { echo "No file $file_list" ; exit 1 ; }
[ "$catdir" ] || { echo "No catalog dir specified" ; exit 2 ; }

go_iter=0
function gonogo() {
    ((go_iter++))
    [ "$((go_iter%NBATCH))" -eq 0 ] && wait
}

#
# starting work
#

ls f??_stack_${cam}.fits > stack_$file_list
[ -s stack_$file_list ] || { echo "No stack files found, exiting..." ; exit 4 ; }

#
# find stars not present in USNO, list in nfile = stack_${cam}_radec.txt
#
run_sex.sh stack_${cam}.fits $inst srcfile=usno_radec.fits -ASSOCSELEC_TYPE -MATCHED & pid0=$!

#
# do photometry on individual files
#
ls `sed -e 's/\.fits/\.resamp\.fits/g' $file_list` 2>/dev/null > r$file_list
x0=`gethead CRPIX1 stack_${cam}.fits`; y0=`gethead CRPIX2 stack_${cam}.fits`
radec2xy.py stack_${cam}.fits $match_file | awk '{print $1-('$x0'),$2-('$y0'),$3,$4,$5}' > match_file.xy

new_files=
for file in `cat r$file_list`; do
    base=${file%'.fits'}
    if [ -d "${base}_dir" ]; then
       echo "Already ran sextractor on $file, skipping..."
    else
        mkdir ${base}_dir
        x=`gethead CRPIX1 $file`; y=`gethead CRPIX2 $file`
        awk '{print $1+'$x',$2+'$y',$3,$4,$5}' match_file.xy > ${base}_dir/sky.list
        new_files="$new_files $file"
    fi
done

for file in $new_files; do
    run_sex.sh $file $inst -ASSOCSELEC_TYPE ALL &
    gonogo
done
wait


#
# find bright stars present in USNO, list in sfile = stack_${cam}_radec0.txt
#
echo find_stars.sh stack_$file_list inst=$inst srcfile=match_file.xy thresh=1.0 cal_mag_max=16.0
find_stars.sh stack_$file_list inst=$inst srcfile=match_file.xy thresh=1.0 cal_mag_max=16.0

# record some statistics
ddoti_getradecrange.py stack_$file_list | sh
for file in `cat stack_$file_list`; do
    base=${file%'.fits'}
    lms=`egrep 'FWHM|10-sigma limiting|Zero' ${base}_dir/${base}.report | awk '{printf("%f ",$4)}END{printf("\n")}' | awk '{lm=$3;if($1<7.3) lm-=log($1/7.3); printf("FWHM=%.2f MAGZERO=%.2f MAGLIM=%.2f\n",$1,$2,lm)}'`
    sethead $lms $file
done

# merge the stellar photometry
sfile=stack_${cam}_radec0.txt
grep '#' f00_stack_${cam}_dir/f00_stack_${cam}_radec.txt > $sfile
for stack_file in `ls f??_stack_${cam}.fits`; do
    x=`gethead CRPIX1 $stack_file`; y=`gethead CRPIX2 $stack_file`
    base=${stack_file%'.fits'}
    grep -v '#' ${base}_dir/${base}_radec.txt.photometry.txt | awk '{print $1,$2,$3,$4,$5,$6,$7,$8-('$x')+'$x0',$9-('$y')+'$y0',$10,$11,$12,$13,$14,$15}'
done >> $sfile

# now back to finding stars not present in USNO, cut on SNR and FWHM
wait $pid0 2>/dev/null
nfile=stack_${cam}_radec.txt
grep -v '#' stack_${cam}_dir/stack_${cam}_radec.txt | awk '{if ($4*'$search_snr'<1) print}' > $nfile
remove_high_props.py $nfile usno_radec.pm.fits

# calibrate the individual frames
nn=`cat match_file.xy | wc -l`
ref=stack_${cam}_radec0.xy
grep -v '#' $sfile | awk '{print $8,$9,$NF,$3,$4}' > $ref

for file in `cat r$file_list`; do
    base=${file%'.fits'}
    calibrate.py ${base}_dir/${base}_radec.txt $ref ${base}_dir/${base}_radec.txt.photometry.txt ${base}_dir/${base}_radec.txt.match.txt 16.0 0.0 > ${base}_dir/${base}.report &
    gonogo
done
wait

# new sources should have fwhm values similar to stars
for stack_file in `ls f??_stack_${cam}.fits`; do
    fwhm=`gethead FWHM $stack_file`
    x=`gethead CRPIX1 $stack_file`; y=`gethead CRPIX2 $stack_file`
    nx=`gethead NAXIS1 $stack_file`; ny=`gethead NAXIS2 $stack_file`
    base=${stack_file%'.fits'}
    mag0=`grep 'Magnitude Offset' ${base}_dir/${base}.report | awk '{print $3;exit}'`
    [ "$mag0" ] || mag0=0.0
    awk '{f=$7; x=$8-('$x0')+'$x'; y=$9-('$y0')+'$y'; if (x>=1 && x<='$nx' && y>=1 && y<='$ny' && f>'$fwhm'/2. && f<'$fwhm'*4) print $1,$2,$3+('$mag0'),$4,$5+('$mag0'),$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' $nfile
done > tmp$nfile
mv tmp$nfile $nfile

# trim out single-epoch sources and clustered detections, limit to ndetect_max
warning=`ddoti_trimcat.py $nfile $match_file $ndetect_max $search_snr`
echo $warning
echo $warning > warning.txt

# now track the new stars in the individual files
echo "# id mag dmag fwhm expos image part" > new_phot_${cam}.txt
grep -v '#' $nfile | awk '{printf("%6d %12.6f %12.6f %8.2f %12.6f %2d 00\n",$NF,$3,$4,$7,$14,0)}' >> new_phot_${cam}.txt
echo "# id ra dec mag dmag fwhm expos image" > cal_phot_${cam}.txt
grep -v '#' $sfile | awk '{printf("%6d %12.8f %12.8f %12.6f %12.6f %8.2f %12.6f %2d\n",$NF,$1,$2,$3,$4,$7,$14,0)}' > cal_phot_${cam}.tmp

[ `cat $file_list | wc -l` -lt 2 ] && do_lightcurves=no
if [ "$do_lightcurves" = "yes" ]; then

    for file in `cat r$file_list`; do
        base=${file%'.fits'}
        base0=${file%'.resamp.fits'}
        n=`grep -n $base0 $file_list | awk -F: '{print $1;exit}'`
        x=`gethead CRPIX1 $file`; y=`gethead CRPIX2 $file`
        nx=`gethead NAXIS1 $file`; ny=`gethead NAXIS2 $file`
        grep -v '#' $nfile | awk '{x=$8-('$x0')+'$x';y=$9-('$y0')+'$y'; if(x>=1 && x<='$nx' && y>=1 && y<='$ny') print x":"y":"$15}' > search_${base}_radec.txt

        nlines=`cat search_${base}_radec.txt | wc -l`
        if [ "$nlines" -gt 0 ]; then
            grep -v '#' ${base}_dir/${base}_radec.txt.photometry.txt | awk '{if($NF==0) print}' > new_${base}_radec.txt
            for src in `cat search_${base}_radec.txt`; do
                x=`echo $src | awk -F: '{print $1;exit}'`
                y=`echo $src | awk -F: '{print $2;exit}'`
                id0=`echo $src | awk -F: '{print $3;exit}'`
                awk 'BEGIN{dis2m=1.e99}{dx=$8-('$x');dy=$9-('$y'); dis2=dx*dx+dy*dy; if(dis2<dis2m) {dis2m=dis2;m=$3;dm=$4;f=$7;dt=$14}}END{if(dis2m<25) printf("%6d %12.6f %12.6f %8.2f %12.6f %2d 00\n",'$id0',m,dm,f,dt,'$n')}' new_${base}_radec.txt >> new_phot_${cam}.txt
            done
        fi
        grep -v '#' ${base}_dir/${base}_radec.txt.photometry.txt | awk '{if($NF>0) {printf("%6d %12.8f %12.8f %12.6f %12.6f %8.2f %12.6f %2d\n",$NF,$1,$2,$3,$4,$7,$14,'$n')}}' >> cal_phot_${cam}.tmp
    done

fi

sort -n -k 1 cal_phot_${cam}.tmp >> cal_phot_${cam}.txt

rm -r f??_stack_${cam}_dir stack_${cam}_dir r$file_list cal_phot_${cam}.tmp 2>/dev/null
