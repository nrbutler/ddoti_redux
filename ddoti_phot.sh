#!/bin/bash

file_list=$1
shift

catdir=
cam=C0

do_psf_plot=no

# source detection
do_lightcurves=yes
ndetect_max=50
search_snr=10

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

if [ -f "${catdir}/usno_radec_${cam}.txt" ]; then
   echo "Using ${catdir}/usno_radec_${cam}.txt"
else
    file0=`head -1 $file_list`
    file0f=`readlink -f $file0`
    here=`pwd`
    cd $catdir
    echo grab_usno_local.sh $file0f
    grab_usno_local.sh $file0f
    cd $here
fi
ln -s ${catdir}/usno_radec_${cam}.txt usno_radec.txt 2>/dev/null
ln -s ${catdir}/usno_radec16_${cam}.txt usno_radec16.txt 2>/dev/null

#
# starting work
#

ls f??_stack_${cam}.fits > stack_$file_list
[ -s stack_$file_list ] || { echo "No stack files found, exiting..." ; exit 4 ; }
ddoti_getradecrange.py stack_$file_list | sh

#
# find stars not present in USNO
#
run_sex.sh stack_${cam}.fits $inst srcfile=usno_radec.txt -ASSOCSELEC_TYPE -MATCHED &
pid0=$!

#
# find bright stars present in USNO
#
pargs=
[ "$do_psf_plot" = "yes" ] && pargs="dao=daoa psf_only=yes"
find_stars.sh stack_$file_list inst=$inst $pargs srcfile=usno_radec16.txt thresh=$match_thresh

for file in `cat stack_$file_list`; do
    base=${file%'.fits'}
    lms=`egrep 'FWHM|10-sigma limiting|Zero' ${base}_dir/${base}.report | awk '{printf("%f ",$4)}END{printf("\n")}' | awk '{lm=$3;if($1<7.3) lm-=log($1/7.3); printf("FWHM=%.2f MAGZERO=%.2f MAGLIM=%.2f\n",$1,$2,lm)}'`
    sethead $lms $file
done
[ "$do_psf_plot" = "yes" ] && psf_plot.py f\*stack_${cam}_dir/psf\*.fits $cam

rm stack_${cam}_radec0.txt 2>/dev/null
x0=`gethead CRPIX1 stack_${cam}.fits`; y0=`gethead CRPIX2 stack_${cam}.fits`
for stack_file in `ls f??_stack_${cam}.fits`; do
    x=`gethead CRPIX1 $stack_file`; y=`gethead CRPIX2 $stack_file`
    base=${file%'.fits'}
    sfile0=${base}_dir/${base}_radec.txt
    grep -v '#' $sfile0 | awk '{print $1,$2,$3,$4,$5,$6,$7,$8-('$x')+'$x0',$9-('$y')+'$y0',$10,$11,$12,$13,$14,$15}' >> stack_${cam}_radec0.txt
done

# now back to finding stars not present in USNO
wait $pid0
sfile=stack_${cam}_dir/stack_${cam}_radec.txt
sfile0=stack_${cam}_dir/stack_${cam}_radec0.txt
grep '#' $sfile > $sfile0
grep -v '#' $sfile | awk '{if ($4*'$search_snr'<1) print}' >> $sfile0
mv $sfile0 $sfile

# new sources should have fwhm values similar to stars
for stack_file in `ls f??_stack_${cam}.fits`; do
    fwhm=`gethead FWHM $stack_file`
    x=`gethead CRPIX1 $stack_file`; y=`gethead CRPIX2 $stack_file`
    nx=`gethead NAXIS1 $stack_file`; ny=`gethead NAXIS2 $stack_file`
    grep -v '#' $sfile | awk '{f=$7; x=$8-('$x0')+'$x'; y=$9-('$y0')+'$y'; if (x>=1 && x<='$nx' && y>=1 && y<='$ny' && f>'$fwhm'/2. && f<'$fwhm'*4) print}'
done > stack_${cam}_radec1.txt

warning=`ddoti_trimcat.py stack_${cam}_radec1.txt $ndetect_max $search_snr`
echo $warning
ndetect=`grep -v '#' stack_${cam}_radec1.txt | wc -l`
grep -v '#' stack_${cam}_radec1.txt | awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,-$15}' >> stack_${cam}_radec0.txt
calibrate.py stack_${cam}_radec0.txt usno_radec16.txt stack_${cam}_radec1.txt stack_${cam}_radec.txt.match.txt > /dev/null
grep '#' stack_${cam}_radec1.txt > stack_${cam}_radec0.txt
grep '#' stack_${cam}_radec1.txt > stack_${cam}_radec.txt
grep -v '#' stack_${cam}_radec1.txt | awk '{if($15>0) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' >> stack_${cam}_radec0.txt
grep -v '#' stack_${cam}_radec1.txt | awk '{if($15<0) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,-$15}' >> stack_${cam}_radec.txt

# now track the new stars in the individual files
echo "# id mag dmag fwhm expos image part" > all_phot_${cam}.txt
awk '{printf("%6d %12.6f %12.6f %8.2f %12.6f %2d 00\n",$NF,$3,$4,$7,$14,0)}' stack_${cam}_radec.txt >> all_phot_${cam}.txt

[ `cat $file_list | wc -l` -lt 2 ] && do_lightcurves=no
if [ "$ndetect" -gt 0 -a "$do_lightcurves" = "yes" ]; then

    ref=stack_${cam}_radec0.txt
    function dophot() {
        local file=$1
        local base=${file%'.fits'}
        local ofile=${file%'.resamp.fits'}_radec.txt
        echo run_sex.sh $file $inst
        run_sex.sh $file $inst > /dev/null 2>&1
        calibrate.py ${base}_dir/${base}_radec.txt $ref > /dev/null
        grep -v '#' ${base}_dir/${base}_radec.txt.photometry.txt | awk '{if($15<0) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,-$15}' > $ofile
    }

    sed -e 's/\.fits/\.resamp\.fits/g' $file_list > r$file_list

    new=stack_${cam}_radec.txt
    x0=`gethead CRPIX1 stack_${cam}.fits`; y0=`gethead CRPIX2 stack_${cam}.fits`
    for file in `cat r$file_list`; do
        base=${file%'.fits'}
        x=`gethead CRPIX1 $file`; y=`gethead CRPIX2 $file`
        nx=`gethead NAXIS1 $file`; ny=`gethead NAXIS2 $file`
        grep -v '#' $new | awk '{x=$8-('$x0')+'$x';y=$9-('$y0')+'$y'; if(x>=1 && x<='$nx' && y>=1 && y<='$ny') print x,y,-$15,$3}' > new_${base}_radec.txt
        n=`cat new_${base}_radec.txt | wc -l`
        if [ "$n" -gt 0 ]; then
            [ -d ${base}_dir ] || mkdir ${base}_dir
            mv new_${base}_radec.txt ${base}_dir/sky.list
            grep -v '#' $ref | awk '{x=$8-('$x0')+'$x';y=$9-('$y0')+'$y'; if(x>=1 && x<='$nx' && y>=1 && y<='$ny') print x,y,NR,$3}' >> ${base}_dir/sky.list
            dophot $file &
            gonogo
        else
            rm new_${base}_radec.txt
        fi
    done
    wait

    n=1
    for file in `cat $file_list`; do
        rfile=${file%'.fits'}_radec.txt
        [ -f "$rfile" ] && awk '{printf("%6d %12.6f %12.6f %8.2f %12.6f %2d 00\n",$NF,$3,$4,$7,$14,'$n')}' $rfile >> all_phot_${cam}.txt
        n=`expr $n + 1`
    done
    rm r$file_list

fi

rm -r f??_stack_${cam}_dir stack_${cam}_dir *${cam}*resamp_dir 2>/dev/null
