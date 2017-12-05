#!/bin/bash

file_list=$1
shift

savedir=`pwd`

# source detection
ndetect_max=50
search_snr=10

# summary
thumb_size=301
thumb_scale=2

t00=`date +%s`

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ "$NBATCH" ] || NBATCH=4

[ -s $file_list ] || { echo "No file $file_list" ; exit 1 ; }

go_iter=0
function gonogo() {
    ((go_iter++))
    [ "$((go_iter%NBATCH))" -eq 0 ] && wait
}

file0=`head -1 $file_list | sed -e 's/\.bz2/\.txt/g'`
base0=${file0%'.fits'}
[ -f $file0 ] || { echo "Cannot find first file $file0" ; exit 2 ; }
cam=`gethead CCD_NAME $file0 | cut -c1-2`
[ "$cam" ] || cam=`basename $file0 | cut -c16-17`
[ "$filter" ] || filter=`gethead FILTER $file0`
[ "$filter" ] || filter=W

dte1=`echo $file0 | sed -e 's/f//g' -e 's/C/ /g' | awk '{print $1}'`
file1=`tail -1 $file_list | sed -e 's/\.bz2/\.txt/g'`
dte2=`echo $file1 | sed -e 's/f//g' -e 's/C/ /g' | awk '{print $1}'`

tag=${cam}_$filter
grb_RA=`gethead CRVAL1 $file0`
grb_DEC=`gethead CRVAL2 $file0`

trig_time=`gethead SEXALT $file0 | sed -e 's/-//g' -e 's/://g'`
[ "$trig_time" ] || trig_time=`gethead DATE-OBS $file0 | sed -e 's/-//g' -e 's/://g'`
gethead -a DATE-OBS EXPTIME @$file_list | awk '{print $2,$3}' | sed -e 's/-//g' -e 's/://g' > times_$file_list
echo ddoti_lc_plots.py all_phot_${cam}.txt stack_${cam}_radec.txt times_$file_list $trig_time
ddoti_lc_plots.py all_phot_${cam}.txt stack_${cam}_radec.txt times_$file_list $trig_time

[ -f stack_${cam}.jpg ] && rm stack_${cam}.jpg
make_circles.sh stack_${cam}.fits stack_${cam}_radec.txt &

#
# make postage stamps
#
file=stack_${cam}.fits
echo "file=$file" > thumb.batch
grep -v '#' stack_${cam}_radec.txt | awk '{x=$8;y=$9;id=$NF;printf("getfits $file %.0f-%.0f %.0f-%.0f -o thumb_%d.fits\n",x-150,x+150,y-150,y+150,id)}' >> thumb.batch
cat thumb.batch | sh >/dev/null 2>&1

x0=`gethead CRPIX1 $file`; y0=`gethead CRPIX2 $file`
dx0=`gethead NAXIS1 $file`; dy0=`gethead NAXIS2 $file`
radec2xy.py stack_${cam}.fits usno_radec.txt > usno_xy.tmp
paste usno_xy.tmp usno_radec.txt | awk '{if($1>=1 && $1<='$dx0' && $2>=1 && $2<='$dy0') printf("%f %f %f_%f\n",$1,$2,$3,$4)}' > usno_xy.txt
rm usno_xy.tmp
n=0
for file in `ls thumb_*fits`; do
    make_circles.sh $file stack_${cam}_radec.txt sz=$thumb_size x0=$x0 y0=$y0 scl=$thumb_scale calfile=usno_xy.txt &
    gonogo
done
wait

sz0=`echo $thumb_size $thumb_scale | awk '{print $1*$2}'`
C11=`gethead CD1_1 stack_${cam}.fits`; C12=`gethead CD1_2 stack_${cam}.fits`
pscale=`echo $C11 $C12 | awk '{print sqrt($1*$1+$2*$2)*3600.}'`
phys_size=`echo $thumb_size $pscale | awk '{printf("%.1f\n",$1*2/60.)}'`
for file in `ls thumb_*_circles.jpg`; do
    file1=`echo $file | sed -e 's/_circles//g'`
    base=${file%'.jpg'}
    base1=${file1%'.jpg'}
    echo "<HTML><BODY><IMG SRC=\"$file\" USEMAP=\"#${base}.map\" WIDTH=\"$sz0\">" > ${base}.html
    echo "<HTML><BODY><IMG SRC=\"$file1\" USEMAP=\"#${base1}.map\" WIDTH=\"$sz0\">" > ${base1}.html
    cat ${base}.map >> ${base}.html
    cat ${base1}.map >> ${base1}.html
    echo "<BR>Physical Size: $phys_size arcmin, Click Red Circle for Lightcurve, Click Elsewsavedir to Toggle Circles<BR></BODY></HTML>" >> ${base}.html
    echo "<BR>Physical Size: $phys_size arcmin, Click to Toggle Circles<BR></BODY></HTML>" >> ${base1}.html
done

#
# compress the stacks
#
rm stack_${cam}.fits.fz stack_${cam}.wt.fits.fz 2>/dev/null
fpack stack_${cam}.fits &
fpack stack_${cam}.wt.fits &
wait

pdir0=photometry_${dte1}_${tag}_${nfiles}
pdir=${savedir}/$pdir0
if [ -d $pdir ]; then
    rm -r ${pdir}/*
else
    mkdir $pdir
fi

#
#make a webpage
#
echo "<meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\">
<HTML><HEAD><TITLE>$dte1 $nfiles $tag redux</TITLE></HEAD>
<BODY BGCOLOR=\"#FFFFFF\" TEXT=\"#003300\">
<FONT SIZE=\"+2\" COLOR=\"#006600\">DDOTI $cam Sources for : RA $grb_RA , Dec $grb_DEC [$nfiles Frame, $dte1 - $dte2]</FONT><P>
<IMG SRC=\"${pdir0}/stack_${dte1}_${tag}_circles.jpg\" USEMAP=\"#stack_${cam}_circles.map\">" > index.html
sed \$d stack_${cam}_circles.map | sed -e "s/lc_/${pdir0}\/thumb_/g" -e "s/\.jpg/_circles\.html/g" >> index.html
echo "</map><BR><A HREF=\"${pdir0}/stack_${dte1}_${tag}.jpg\" target=\"_blank\">Image without circles</A>" >> index.html
echo "<P><FONT SIZE=\"+2\" COLOR=\"#006600\">USNO-B1 R Band Photometry (Uncatalogued Sources, SNR>$search_snr, Nmax=$ndetect_max):</FONT><BR><P><PRE>" >> index.html
grep '#' all_phot_${cam}.txt.summary.txt >> index.html
grep -v '#' all_phot_${cam}.txt.summary.txt | eval awk "'{printf(\"<A HREF=\\\"$pdir0/thumb_%s_circles.html\\\" target=\\\"_blank\\\">%s</A>\n\",\$1,\$1)}'" > phot$$.tmp0
grep -v '#' all_phot_${cam}.txt.summary.txt | cut -c8- > phot$$.tmp1
paste phot$$.tmp0 phot$$.tmp1 >> index.html
echo "</PRE><BR>${warning}<HR>Substacks made for this field:<BR><PRE>
# File              RA1          RA2          DEC1         DEC2        EXPOSURE  FWHM   Mag-10s  ZeroPoint" >> index.html
gethead -a RA0 RA1 DEC0 DEC1 EXPTIME FWHM MAGLIM MAGZERO @stack_$file_list | awk '{printf("%s %12.8f %12.8f %12.8f %12.8f %8.1f %6.1f %8.2f %8.2f\n",$1,$2,$3,$4,$5,$6,$7,$8,$9)}' >> index.html
echo "</PRE>" >> index.html
[ -f "${cam}_psf_plot.jpg" ] && echo "<A HREF=\"${pdir0}/stack_${dte1}_${tag}_psf_plot.jpg\" target=\"_blank\">PSF Summary</A>" >> index.html
[ -f "${base0}_dist.jpg" ] && echo "<A HREF=\"${pdir0}/stack_${dte1}_${tag}_dist.jpg\" target=\"_blank\">Distortion Correction</A>" >> index.html

task_tot=`date +%s | awk '{printf("Total Execution Time for All Tasks: %.2f minutes\n",($1-'$t00')/60.)}'`
upt=`uptime`
echo $task_tot
echo $upt

echo "<HR WIDTH=\"100%\"> $task_tot <A HREF=\"stack_${dte1}_${tag}_${nfiles}.log\" target=\"_blank\">logfile</A><BR>
$upt <BR> Last Updated: $dte (natbutler@asu.edu)</BODY></HTML>" >> index.html


#
# ship out all the summary files
#
mv index.html ${savedir}/stack_${dte1}_${tag}_${nfiles}.html
rm ${savedir}/current.html 2>/dev/null
ln -s ${savedir}/stack_${dte1}_${tag}_${nfiles}.html ${savedir}/current.html

mv all_phot_${cam}.txt ${pdir}/all_phot_${dte1}_${tag}.txt
mv stack_${cam}_radec.txt ${pdir}/stack_${dte1}_${tag}_radec.txt
[ -f "${cam}_psf_plot.jpg" ] && cp ${cam}_psf_plot.jpg ${pdir}/stack_${dte1}_${tag}_psf_plot.jpg
cp ${base0}_dist.jpg ${pdir}/stack_${dte1}_${tag}_dist.jpg
mv all_phot_${cam}.txt.summary.txt ${pdir}/all_phot_${dte1}_${tag}_summary.txt
mv lc_*.jpg thumb_*.jpg thumb_*.html $pdir 2>/dev/null
mv stack_${cam}.jpg ${pdir}/stack_${dte1}_${tag}.jpg
mv stack_${cam}_circles.jpg ${pdir}/stack_${dte1}_${tag}_circles.jpg
mv stack_${cam}.fits.fz ${savedir}/stack_${dte1}_${tag}.fits.fz
mv stack_${cam}.wt.fits.fz ${savedir}/stack_${dte1}_${tag}.wt.fits.fz

rm thumb.batch phot$$.tmp0 phot$$.tmp1 thumb_*.fits thumb_*.map stack_*.map times_$file_list usno_xy.txt 2>/dev/null
