#!/bin/sh

fitsfile=$1
shift
radeclist=$1
shift

sz=1024
scl=1

x0=`gethead CRPIX1 $fitsfile`
y0=`gethead CRPIX2 $fitsfile`
calfile=

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

x=`gethead CRPIX1 $fitsfile`
y=`gethead CRPIX2 $fitsfile`

nx=`gethead NAXIS1 $fitsfile`
ny=`gethead NAXIS2 $fitsfile`
sx=$sz
sy=`echo "$sz $nx $ny" | awk '{printf("%.0f\n",$1*$3/$2)}'`
#sx=`echo "$sz $nx $ny" | awk '{if ($2>$3) print $1; else printf("%.0f\n",$1*$2/$3)}'`
#sy=`echo "$sz $nx $ny" | awk '{if ($3>$2) print $1; else printf("%.0f\n",$1*$3/$2)}'`

facx=`echo $sx $nx | awk '{printf("%.6f\n",$2/$1)}'`
facy=`echo $sy $ny | awk '{printf("%.6f\n",$2/$1)}'`

base=${fitsfile%'.fits'}
[ -f ${base}.jpg ] || fits2jpg.py $fitsfile $sz

echo "<map name=\"${base}_circles.map\">" > ${base}_circles.map
echo "<map name=\"${base}.map\">" > ${base}.map

redcircles=
if [ -f "$radeclist" ]; then
    grep -v '#' $radeclist | awk '{x=$8-('$x0')+'$x';y=$9-('$y0')+'$y'; if (x>=1 && x<='$nx' && y>=1 && y<='$ny') print x,y,$3,$NF}' > radec$$.txt

    redcircles=`awk '{x=($1-1)/'$facx';y=('$ny'-$2)/'$facy';n=$NF; printf(" -stroke red -fill None -draw '\''circle %.0f,%.0f %.0f,%0.f'\'' -stroke black -fill yellow -draw '\''translate %.0f,%.0f text -5,-11 \"%d\"'\''",x,y,x,y+5,x,y,n)}' radec$$.txt`

    awk '{x=($1-1)/'$facx';y=('$ny'-$2)/'$facy';id=$NF; printf("<area shape=\"circle\" coords=\"%.0f,%.0f,10\" href=\"lc_%d.jpg\" target=\"_blank\" title=\"Source %d, R=%.1f\">\n",x*'$scl',y*'$scl',id,id,$3)}' radec$$.txt > ${base}.tmp
    cat ${base}.tmp >> ${base}_circles.map
    cat ${base}.tmp >> ${base}.map

    rm radec$$.txt ${base}.tmp
fi

bluecircles=
if [ "$calfile" ]; then
    if [ -f "$calfile" ]; then
        grep -v '#' $calfile | awk '{x=$1-('$x0')+'$x';y=$2-('$y0')+'$y'; if (x>=1 && x<='$nx' && y>=1 && y<='$ny') print x,y,$NF}' > radec$$.txt
        bluecircles=`awk '{x=($1-1)/'$facx';y=('$ny'-$2)/'$facy'; printf(" -stroke blue -fill None -draw '\''circle %.0f,%.0f %.0f,%0.f'\'' ",x,y,x,y+2)}' radec$$.txt`
        awk '{x=($1-1)/'$facx';y=('$ny'-$2)/'$facy'; printf("<area shape=\"circle\" coords=\"%.0f,%.0f,10\" target=\"_blank\" title=\"USNO %s\">\n",x*'$scl',y*'$scl',$3)}' radec$$.txt > ${base}.tmp
        cat ${base}.tmp >> ${base}_circles.map
        cat ${base}.tmp >> ${base}.map
        rm radec$$.txt ${base}.tmp
    fi
fi

sz0=`echo $sz $scl | awk '{print $1*$2}'`
echo "<area shape=\"rect\" coords=\"0,0,$sz0,$sz0\" href=\"${base}_circles.html\" title=\"Click to Return Circles\"></map>" >> ${base}.map
echo "<area shape=\"rect\" coords=\"0,0,$sz0,$sz0\" href=\"${base}.html\" title=\"Click to Remove Circles\"></map>" >> ${base}_circles.map

cmd="convert -quality 75 ${base}.jpg -pointsize 16 $redcircles $bluecircles ${base}_circles.jpg"
eval "$cmd"
