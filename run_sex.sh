#!/bin/sh

file=`readlink -f $1`
shift

inst=$1
if [ "$inst" ]; then
   shift
else
    inst=nircam
fi

do_dao=no
if [ "$1" = "dao" ]; then
    do_dao=yes
    shift
fi
aperture=
if [ "$1" = "daoa" ]; then
    aperture="--aperture"
    do_dao=yes
    shift
fi
if [ "$1" = "daoan" ]; then
    aperture="--aperture --noremove_contam"
    do_dao=yes
    shift
fi
dual=no
if [ "$1" = "dual" ]; then
    dual=yes
    shift
fi
psf_only=no
psfonly=
if [ "$1" = "psfonly" ]; then
    psf_only=yes
    psfonly="--psfonly"
    shift
fi

srcfile=`echo $1 | awk -F= '/srcfile/{print $2}'`
[ "$srcfile" ] && shift

psffile=`echo $1 | awk -F= '/psffile/{print $2}'`
[ "$psffile" ] && shift

extras=$@

rms_file=`echo $file | sed -e 's/\.fits/\.rms\.fits/g'`
[ -f "$rms_file" ] || {
    wt_file=`echo $file | sed -e 's/\.fits/\.wt\.fits/g'`
    [ -f "$wt_file" ] || wt_file=`echo $file | sed -e 's/\.fits/\.weight\.fits/g'`
}
flag_file=`echo $file | sed -e 's/\.fits/\.flag\.fits/g'`

file0=`basename $file`
base=${file0%'.fits'}
base_dir=`dirname $file`
dir=${base}_dir
[ -d $dir ] || mkdir $dir

[ "$psffile" ] && cp $psffile $dir
cd $dir

[ "$inst" ] || inst=default
sx=`gethead NAXIS1 $file`
sy=`gethead NAXIS2 $file`

do_assoc=no
if [ "$srcfile" ]; then
    [ -f "${base_dir}/$srcfile" ] && grep -v '#' ${base_dir}/$srcfile | awk '{print $1,$2,NR,$3}' > sky_radec.list
fi
if [ -f sky.list -o -f sky_radec.list ]; then
    echo "Running in assoc mode using sky.list file"
    do_assoc=yes
    inst=${inst}_assoc
    if [ -f sky_radec.list ]; then
        radec2xy.py $file sky_radec.list > sky_xy.list
        awk '{print $3,$4}' sky_radec.list > sky_xy.tmp
        paste sky_xy.list sky_xy.tmp | awk '{if($1>=1 && $1<='$sx' && $2>=1 && $2<='$sy') print}' > sky.list
        rm sky_xy.list sky_xy.tmp
        ns=`cat sky.list | wc -l`
        [ "$ns" -le 0 ] && { echo "No sources present to match" ; exit 1 ; }
    fi
fi

[ -f "$flag_file" ] && ln -s $flag_file flag.fits

gain=`gethead GAIN $file`
[ "$gain" ] || {
    gain=`gethead EGAIN $file`
    [ "$gain" ] && sethead GAIN=$gain $file
    [ "$gain" ] || gain=1.0
}

wfile=
if [ -f "$rms_file" ];then
    extras="-WEIGHT_TYPE MAP_RMS,MAP_RMS -WEIGHT_GAIN N,N -WEIGHT_IMAGE $rms_file,$rms_file -GAIN 1.e9 $extras"
    wfile=$rms_file
elif [ -f "$wt_file" ];then
    extras="-WEIGHT_TYPE MAP_WEIGHT,MAP_WEIGHT -WEIGHT_GAIN Y,Y -WEIGHT_IMAGE $wt_file,$wt_file -WEIGHT_THRESH 0 -GAIN $gain $extras"
    wfile=$wt_file
else
    extras="-WEIGHT_TYPE NONE,NONE $extras"
fi

if [ "$do_dao" = "yes" ]; then
    extras="$extras -CHECKIMAGE_TYPE OBJECTS -CHECKIMAGE_NAME mask.fits"
fi

bs=`gethead BACK_SZ $file`
[ "$bs" ] && extras="$extras -BACK_SIZE $bs"

t0=`gethead DATE-OBS $file | sed -e 's/-//g' -e 's/://g' -e 's/T/_/g'`
t1=`gethead DATE-OBE $file | sed -e 's/-//g' -e 's/://g' -e 's/T/_/g'`
[ "$t1" ] || t1=$t0

ccd=`gethead CCD_NAME $file`
sat_level=`gethead SATURATE $file`
[ "$sat_level" ] || sat_level=60000.0
extras="-SATUR_LEVEL $sat_level $extras -CATALOG_NAME test.cat "

echo sextractor -c ${SEXTRACTOR_DIR}/${inst}.sex $file $extras -PARAMETERS_NAME ${SEXTRACTOR_DIR}/${inst}.param -FILTER_NAME ${SEXTRACTOR_DIR}/${inst}.conv -STARNNW_NAME ${SEXTRACTOR_DIR}/${inst}.nnw
sextractor -c ${SEXTRACTOR_DIR}/${inst}.sex $file $extras -PARAMETERS_NAME ${SEXTRACTOR_DIR}/${inst}.param -FILTER_NAME ${SEXTRACTOR_DIR}/${inst}.conv -STARNNW_NAME ${SEXTRACTOR_DIR}/${inst}.nnw

exposure=`gethead EXPTIME0 $file`
[ "$exposure" ] || exposure=`gethead EXPTIME $file`
am=`gethead AIRMASS $file`
[ "$am" ] && airmass=`gethead AM_COEF $file | awk '{if('$am'>1) print $1*('$am'-1); else print 0.0}'`
[ "$airmass" ] || airmass=0.0

nsources=`grep -v '#' test.cat | wc -l`
if [ "$nsources" -eq 0 ]; then
    echo "Found no sources..."
    echo "# RA DEC mag dmag mag_big dmag_big FWHM x y xa ya x2a y2a expos num (ccd $ccd gain $gain exposure $exposure airmass_ext $airmass sex_zero 25.0 t0 $t0 t1 $t1 )" > ${base}_radec.txt
fi

[ -f test_assoc.cat ] && mv test_assoc.cat test.cat

echo "# RA DEC mag dmag mag_big dmag_big FWHM x y xa ya x2a y2a expos num (ccd $ccd gain $gain exposure $exposure airmass_ext $airmass sex_zero 25.0 t0 $t0 t1 $t1 )" > ${base}_radec.txt
n=`grep -v '#' test.cat | wc -l`
if [ "$n" -gt 0 ]; then
    grep -v '#' test.cat | awk '{print $1,$2,$3}' > xy_in.txt
    grep -v '#' test.cat | awk '{if (NF==13) print $11,$12,$4,$6,$5,$7,$13,$1,$2; else print $25,$26,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$27,$1,$2}' >> ${base}_radec.txt

    # replace x and y positions with relative positions
    grep '#' ${base}_radec.txt > ${base}_radec.tmp
    grep -v '#' ${base}_radec.txt > ${base}_radec.tmp0
    awk '{x=$1/'$sx';y=$2/'$sy'; print x,y,0.5*(3*x*x-1),0.5*(3*y*y-1),'$exposure',$3}' xy_in.txt > xy_out.txt
    paste -d ' ' ${base}_radec.tmp0 xy_out.txt >> ${base}_radec.tmp
    grep -v 'nan' ${base}_radec.tmp > ${base}_radec.txt
    rm ${base}_radec.tmp0 ${base}_radec.tmp
fi

if [ "$do_dao" = "yes" ]; then
    cp ${base}_radec.txt ${base}_radec0.txt
    if [ -f "${base_dir}/force_phot.txt" ]; then
        awk '{print $1,$2,-NR}' ${base_dir}/force_phot.txt > xy_in.txt
        if [ -s xy_in.txt ]; then
            echo "Forcing photometry at positions in ${base_dir}/force_phot.txt:"
            cat xy_in.txt
            radec2xy.py $file ${base_dir}/force_phot.txt > xy_force.txt
            awk '{x=$1/'$sx';y=$2/'$sy'; print x,y,0.5*(3*x*x-1),0.5*(3*y*y-1),'$exposure',-NR}' xy_force.txt > xy_out.txt
            dat=`grep -v '#' ${base}_radec.txt | awk '{printf("%f,%f,%f,%f,%f\n",$3,$4,$5,$6,$7);exit}'`
            [ "$dat" ] || dat="19.0,0.1,19.0,0.1,10.0"
            paste xy_in.txt xy_force.txt xy_out.txt | eval awk "'{print \$1,\$2,$dat,\$4,\$5,\$6,\$7,\$8,\$9,\$10,\$3}'" >> ${base}_radec.txt
        fi
    fi
    echo sources2mask.py $wfile mask.fits mask${base}.fits $fwhm
    sources2mask.py $wfile mask.fits mask${base}.fits $fwhm
    echo psf_fit.py $file ${base}_radec.txt $psfonly --satlevel $sat_level $aperture > ${base}_radec1.txt
    psf_fit.py $file ${base}_radec.txt $psfonly --satlevel $sat_level $aperture > ${base}_radec1.txt
    [ -f ${base}_psf_stats.txt ] && mv ${base}_psf_stats.txt psf_stats.txt
    if [ -s ${base}_radec1.txt -a "$psf_only" = "no" ]; then
        mv ${base}_radec1.txt ${base}_radec.txt
    else
        echo "PSF fitting not applied or failed."
    fi
fi
