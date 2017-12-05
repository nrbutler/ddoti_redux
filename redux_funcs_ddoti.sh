#!/bin/bash

export REDUX_BASE_DIR=/data_storeA/ddoti
export PATH=${REDUX_BASE_DIR}/bin:/usr/local/astrometry/bin:${REDUX_BASE_DIR}/python_modules:/usr/bin:/usr/local/bin:/bin
export PYTHONPATH=${REDUX_BASE_DIR}/python_modules

export REDUX_LOCKFILE=${REDUX_BASE_DIR}/ddoti.lock
export SEXTRACTOR_DIR=${REDUX_BASE_DIR}/sextractor
export SWARP_DIR=${REDUX_BASE_DIR}/swarp
export ASTNET_DIR=${REDUX_BASE_DIR}/astnet
export CALFILE_DIR=${REDUX_BASE_DIR}/calfiles

# set test=echo for testing purposes
#test=echo
test=

cleanup=yes

# all functions look at a day's worth of data defined by TODAY:
export TODAY=`date -u +20%y%m%d`

# location of raw data and summary (web) files
export data_pc=tcs-a
export raw_data_dir=/images/test
export raw_data_archive=/usr/local/var/archive/date
export raw_data_archive0=/archive0
export web_pc=tcs-a
export web_user=reducer
if [ "$HOSTNAME" = "ratir" -o "$HOSTNAME" = "legion" ]; then
    export data_pc=ratir
    export web_user=nrbutler
    export web_pc=butler.lab.asu.edu
fi
export web_data_dir=/usr/local/var/www/main/ddoti
[ "$HOSTNAME" = "ratir" -o "$HOSTNAME" = "legion" ] && export web_data_dir=public_html/ratir_testing
echo "data pc: ${data_pc}"
echo "web_site pc: ${web_user}@${web_pc}"
export grb_RA=
export grb_DEC=

export NBATCH=`grep processor /proc/cpuinfo | wc -l`

cd $REDUX_BASE_DIR

function ratir_setdirs() {
    # default directories, can be changed manually
    $test rsync -a -f"+ */" -f"- *" --chmod="g=rwx" ${data_pc}:${raw_data_dir}/ $REDUX_BASE_DIR/ 2>/dev/null
    export BIAS_DIRS=`ls -d ${TODAY}/20*-0006/* 2>/dev/null`
    export DARK_DIRS=`ls -d ${TODAY}/20*-0007/* 2>/dev/null`
    export FLAT_DIRS=`ls -d ${TODAY}/20*-0008/* 2>/dev/null`
    export NIR_STANDARD_DIRS=`ls -d ${TODAY}/20*-0000/* 2>/dev/null`
    export OPT_STANDARD_DIRS=`ls -d ${TODAY}/20*-0001/* 2>/dev/null`
    export GRB_DIRS=`ls -d ${TODAY}/2012A-1000/* 2>/dev/null`
    export LIGO_DIRS=`ls -d ${TODAY}/2015B-1000/* ${TODAY}/2017A-1000/* 2>/dev/null`
    export AST_DIRS=`ls -d ${TODAY}/2014B-1004/* ${TODAY}/2015A-1005/* ${TODAY}/2015B-1013/* ${TODAY}/2016A-1012/* ${TODAY}/2016B-1004/* ${TODAY}/2017A-1009/* ${TODAY}/2017B-1005/* 2>/dev/null`
    env | grep _DIRS
    echo "set source_list=<one of these> to activate"
    echo "for a GRB, also set grb_RA and grb_DEC (decimal)"
    env | grep 'grb_'
}

function ratir_copy_files() {
    # populate a data directory with raw images from the data server
    for dir in $source_list; do
        $test rsync -avz --exclude '*g.fits' ${data_pc}:${raw_data_dir}/$dir/20*fits $dir/
    done
}

function ratir_copy_archive_files() {
    # populate a data directory with raw images from the data server
    for dir in $source_list; do
        $test rsync -azv --exclude '*g.fits*' ${data_pc}:${raw_data_archive}/$dir/20*fits.* $dir/
    done
}

function ratir_copy_archive0_files() {
    # populate a data directory with raw images from the data server
    for dir in $source_list; do
        $test rsync -azv --exclude '*g.fits.*' ${data_pc}:${raw_data_archive0}/$dir/20*fits.* $dir/
    done
}

if [ "$HOSTNAME" = "ratir" -o "$HOSTNAME" = "legion" ]; then
    function ratir_copy_files() {
        for dir in $source_list; do
            $test ln -s ${raw_data_archive}/$dir/20*fits.* $dir/ 2>/dev/null
            $test rm ${dir}/*g.fits.* 2>/dev/null
        done
    }
fi

function merge_flat_dirs() {
    ndirs=`echo $FLAT_DIRS | wc -w`
    if [ "$ndirs" -gt 1 ]; then
        dir0=`echo $FLAT_DIRS | awk '{print $1}'`
        for dir in `echo $FLAT_DIRS | awk '{for(i=2;i<=NF;i++) print $i}'`; do
            mv $dir/20*fits* $dir0
            rmdir $dir
        done
        FLAT_DIRS=$dir0
    fi
}

function make_lists() {
    # simple listing function used in manifest creation below
    # we only require the FILTER keyword to be present for C0
    ( gethead -auc FILTER `ls 20*fits.txt 2>/dev/null` 2>/dev/null | sed -e 's/\.txt/\.bz2/g' -e 's/\///g' ;
    gethead -au FILTER `ls 20*fits 2>/dev/null` 2>/dev/null | sed -e 's/\///g' ) | sort -u > manifest$$.txt
    awk '{ccd=substr($1,16,2); f=""; if (ccd~/C0/) f=""$2"_"; print $1>ccd"_"f"list.txt"}' manifest$$.txt
    rm manifest$$.txt C0_____list.txt C0_UNK_list.txt C0_list.txt 2>/dev/null
    cat C0*list.txt 2>/dev/null > C0_list.txt
}

function ratir_create_manifest() {
    # group the raw fits files for reduction
    for dir in $source_list; do
        n=`ls $dir | grep fits | wc -l`
        [ "$n" -eq 0 ] && continue
        $test cd $dir
        $test make_lists
        $test cd $REDUX_BASE_DIR
    done
}

function ratir_redux() {
    # science data reduction (after bias and flat construction), work in parallel
    if [ -f "$REDUX_LOCKFILE" ]; then
        echo "lockfile $REDUX_LOCKFILE present, aborting..."
    else
        $test touch $REDUX_LOCKFILE
        for dir in $source_list; do
            n=`ls $dir | grep fits | wc -l`
            [ "$n" -eq 0 ] && continue
            TODAY=`echo $dir | awk -F/ '{print $1}'`
            cd $dir
            $test get_grb_info.sh logfile=redux.log ra0=$grb_RA dec0=$grb_DEC
            echo "Reducing C0/C1 in parallel..."
            for list in `ls C0_*_list.txt C1_list.txt 2>/dev/null`; do
                $test wrap_redux_h2rg $list ngroups=300 discard_n_lt_ngroups=no logfile=redux.log ra0=$grb_RA dec0=$grb_DEC cleanup=$cleanup $redux_h2rg_args &
            done
            [ "$HOSTNAME" != "legion" ] && wait
            echo "Reducing C2 and C3 in parallel..."
            for list in `ls C[2,3]_list.txt 2>/dev/null`; do
                $test wrap_redux_h2rg $list ngroups=300 discard_n_lt_ngroups=no logfile=redux.log ra0=$grb_RA dec0=$grb_DEC cleanup=$cleanup $redux_h2rg_args &
            done
            wait
            [ -f C4_list.txt ] && $test wrap_redux_h2rg C4_list.txt match_cam=C1 ngroups=300 discard_n_lt_ngroups=no scale_flux=no scale_flux_indv=no logfile=redux.log ra0=$grb_RA dec0=$grb_DEC cleanup=$cleanup $redux_h2rg_args
            echo "Generating summary web-page..."
            touch new_data.txt
            $test get_imstats.sh logfile=redux.log ra0=$grb_RA dec0=$grb_DEC $redux_h2rg_args
            cd $REDUX_BASE_DIR
            $test ratir_copy2web
        done
        $test rm $REDUX_LOCKFILE
    fi
}

function ratir_ligo_redux() {
    # science data reduction (after bias and flat construction), work in parallel
    if [ -f "$REDUX_LOCKFILE" ]; then
        echo "lockfile $REDUX_LOCKFILE present, aborting..."
    else
        $test touch $REDUX_LOCKFILE
        redux_h2rg_args="sky_only=yes split_back=no"
        for dir in $source_list; do
            n=`ls $dir | grep fits | wc -l`
            [ "$n" -eq 0 ] && continue
            TODAY=`echo $dir | awk -F/ '{print $1}'`
            cd $dir
            $test get_grb_info.sh logfile=redux.log ra0=$grb_RA dec0=$grb_DEC
            echo "Reducing C0/C1 in parallel..."
            for list in `ls C0_*_list.txt C1_list.txt 2>/dev/null`; do
                $test wrap_redux_h2rg $list ngroups=300 discard_n_lt_ngroups=no logfile=redux.log ra0=$grb_RA dec0=$grb_DEC cleanup=$cleanup $redux_h2rg_args &
            done
            [ "$HOSTNAME" != "legion" ] && wait
            echo "Reducing C2 and C3 in parallel..."
            for list in `ls C[2,3]_list.txt 2>/dev/null`; do
                $test wrap_redux_h2rg $list ngroups=300 discard_n_lt_ngroups=no logfile=redux.log ra0=$grb_RA dec0=$grb_DEC cleanup=$cleanup $redux_h2rg_args &
            done
            wait
            [ -f C4_list.txt ] && $test wrap_redux_h2rg C4_list.txt match_cam=C1 ngroups=300 discard_n_lt_ngroups=no scale_flux=no scale_flux_indv=no logfile=redux.log ra0=$grb_RA dec0=$grb_DEC cleanup=$cleanup $redux_h2rg_args
            cd $REDUX_BASE_DIR
        done

        rm sky_bank/sky_${TODAY}*.fits 2>/dev/null
        for cam in C0 C1 C2 C3 C4; do
            for dir in $source_list; do
                $test ls ${dir}/v*${cam}*/skystack_${cam}*.fits >> sky_bank/sky_${TODAY}_${cam}.txt
            done
            [ -s sky_bank/sky_${TODAY}_${cam}.txt ] && $test back_combine.py sky_bank/sky_${TODAY}_${cam}.txt sky_bank/sky_${TODAY}_${cam}.fits
            $test rm sky_bank/sky_${TODAY}_${cam}.txt sky_bank/expmapsky_${TODAY}_${cam}.fits 2>/dev/null
        done

        redux_h2rg_args="match_cam=self"
        for dir in $source_list; do
            n=`ls $dir | grep fits | wc -l`
            [ "$n" -eq 0 ] && continue
            cd $dir
            echo "Reducing C0/C1 in parallel..."
            for list in `ls C0_*_list.txt C1_list.txt 2>/dev/null`; do
                cam=`echo $list | cut -c1-2`
                $test wrap_redux_h2rg $list ngroups=300 discard_n_lt_ngroups=no logfile=redux.log ra0=$grb_RA dec0=$grb_DEC cleanup=$cleanup $redux_h2rg_args skyfile=${REDUX_BASE_DIR}/sky_bank/sky_${TODAY}_${cam}.fits regrind=yes &
            done
            [ "$HOSTNAME" != "legion" ] && wait
            echo "Reducing C2 and C3 in parallel..."
            for list in `ls C[2,3]_list.txt 2>/dev/null`; do
                cam=`echo $list | cut -c1-2`
                $test wrap_redux_h2rg $list ngroups=300 discard_n_lt_ngroups=no logfile=redux.log ra0=$grb_RA dec0=$grb_DEC cleanup=$cleanup $redux_h2rg_args skyfile=${REDUX_BASE_DIR}/sky_bank/sky_${TODAY}_${cam}.fits regrind=yes &
            done
            wait
            [ -f C4_list.txt ] && $test wrap_redux_h2rg C4_list.txt match_cam=C1 ngroups=300 discard_n_lt_ngroups=no scale_flux=no scale_flux_indv=no logfile=redux.log ra0=$grb_RA dec0=$grb_DEC cleanup=$cleanup $redux_h2rg_args skyfile=${REDUX_BASE_DIR}/sky_bank/sky_${TODAY}_C4.fits regrind=yes
            echo "Generating summary web-page..."
            touch new_data.txt
            $test get_imstats.sh logfile=redux.log ra0=$grb_RA dec0=$grb_DEC $redux_h2rg_args
            cd $REDUX_BASE_DIR
            $test ratir_copy2web
        done
        redux_h2rg_args=
        $test rm $REDUX_LOCKFILE
    fi
}

function ratir_ast_redux() {
    # science data reduction (after bias and flat construction), work in parallel
    if [ -f "$REDUX_LOCKFILE" ]; then
        echo "lockfile $REDUX_LOCKFILE present, aborting..."
    else
        $test touch $REDUX_LOCKFILE
        for dir in $source_list; do
            n=`ls $dir | grep fits | wc -l`
            [ "$n" -eq 0 ] && continue
            TODAY=`echo $dir | awk -F/ '{print $1}'`
            cd $dir
            $test get_horizons.sh
            $test get_grb_info.sh logfile=redux.log ra0=$grb_RA dec0=$grb_DEC
            echo "Reducing C0/C1 in parallel..."
            for list in `ls C0_*_list.txt C1_list.txt 2>/dev/null`; do
                $test wrap_redux_h2rg $list ngroups=300 discard_n_lt_ngroups=no logfile=redux.log ra0=$grb_RA dec0=$grb_DEC cleanup=$cleanup do_asteroid=yes $redux_h2rg_args &
            done
            [ "$HOSTNAME" != "legion" ] && wait
            echo "Reducing C2 and C3 in parallel..."
            for list in `ls C[2,3]_list.txt 2>/dev/null`; do
                $test wrap_redux_h2rg $list ngroups=300 discard_n_lt_ngroups=no logfile=redux.log ra0=$grb_RA dec0=$grb_DEC cleanup=$cleanup do_asteroid=yes $redux_h2rg_args &
            done
            wait
            [ -f C4_list.txt ] && $test wrap_redux_h2rg C4_list.txt match_cam=C1 ngroups=300 discard_n_lt_ngroups=no scale_flux=no scale_flux_indv=no logfile=redux.log ra0=$grb_RA dec0=$grb_DEC cleanup=$cleanup do_asteroid=yes $redux_h2rg_args
            echo "Generating summary web-page..."
            touch new_data.txt
            $test get_imstats.sh logfile=redux.log ra0=$grb_RA dec0=$grb_DEC $redux_h2rg_args
            cd $REDUX_BASE_DIR
            $test ratir_copy2web
        done
        $test rm $REDUX_LOCKFILE
    fi
}

function ratir_grb_redux() {
    # science data reduction (after bias and flat construction), work in parallel
    if [ -f "$REDUX_LOCKFILE" ]; then
        echo "lockfile $REDUX_LOCKFILE present, aborting..."
    else
        $test touch $REDUX_LOCKFILE
        for dir in $source_list; do
            n=`ls $dir | grep fits | wc -l`
            [ "$n" -eq 0 ] && continue
            TODAY=`echo $dir | awk -F/ '{print $1}'`
            cd $dir
            logfile=grb_redux.log
            [ -f "$logfile" ] && {
                n=`ls grb_redux_*_.log 2>/dev/null | sort -r | awk -F_ '{print $3+1;exit}'`
                [ "$n" ] || n=0
                $test mv $logfile grb_redux_${n}_.log
            }
            $test get_grb_info.sh logfile=$logfile ra0=$grb_RA dec0=$grb_DEC
            echo "Reducing C0/C1 in parallel..."
            for list in `ls C0_*_list.txt C1_list.txt 2>/dev/null`; do
                $test wrap_redux_h2rg $list logfile=$logfile ra0=$grb_RA dec0=$grb_DEC cleanup=$cleanup $redux_h2rg_args &
            done
            if [ "$HOSTNAME" != "legion" ]; then
                wait
                cd $REDUX_BASE_DIR
                $test ratir_copy2web
                cd $dir
                [ -d version0 ] || {
                    echo "Generating summary web-page..."
                    $test get_imstats.sh logfile=$logfile ra0=$grb_RA dec0=$grb_DEC $redux_h2rg_args
                    cd $REDUX_BASE_DIR
                    $test ratir_copy2web
                    cd $dir
                }
            fi
            echo "Reducing C2/C3 in parallel.."
            for list in `ls C[2,3]_list.txt 2>/dev/null`; do
                $test wrap_redux_h2rg $list logfile=$logfile ra0=$grb_RA dec0=$grb_DEC cleanup=$cleanup $redux_h2rg_args &
            done
            wait
            [ -f C4_list.txt ] && $test wrap_redux_h2rg C4_list.txt match_cam=C1 scale_flux=no scale_flux_indv=no logfile=$logfile ra0=$grb_RA dec0=$grb_DEC cleanup=$cleanup $redux_h2rg_args
            echo "Generating summary web-page..."
            $test get_imstats.sh logfile=$logfile ra0=$grb_RA dec0=$grb_DEC $redux_h2rg_args
            cd $REDUX_BASE_DIR
            $test ratir_copy2web
        done
        $test rm $REDUX_LOCKFILE
    fi
}

function ratir_bias() {
    # create bias frames and store to bias bank, work in parallel
    for dir in $source_list; do
        n=`ls $dir | grep fits | wc -l`
        [ "$n" -eq 0 ] && continue
        cd $dir
        for list in `ls C[0,1,4]_list.txt 2>/dev/null`; do
            echo bias_opt $list cleanup=$cleanup
            $test bias_opt $list cleanup=$cleanup &
        done
        wait
        cd $REDUX_BASE_DIR
    done
}

function ratir_dark() {
    # create dark frames and store to dark bank, work in parallel
    for dir in $source_list; do
        n=`ls $dir | grep fits | wc -l`
        [ "$n" -eq 0 ] && continue
        cd $dir
        for list in `ls C[0,1,4]_list.txt 2>/dev/null`; do
            echo dark_opt $list cleanup=$cleanup
            $test dark_opt $list cleanup=$cleanup &
        done
        wait
        cd $REDUX_BASE_DIR
    done
}

function ratir_flat() {
    # create flats frames and store to flat bank, work in parallel
    for dir in $source_list; do
        n=`ls $dir | grep fits | wc -l`
        [ "$n" -eq 0 ] && continue
        cd $dir
        for list in `ls C0_*_list.txt C1_list.txt C4_list.txt 2>/dev/null`; do
            echo flat_h2rg $list cleanup=$cleanup
            $test flat_h2rg $list cleanup=$cleanup &
        done
        wait
        for list in `ls C[2,3]_list.txt 2>/dev/null`; do
            echo $test flat_h2rg $list cleanup=$cleanup
            $test flat_h2rg $list cleanup=$cleanup &
        done
        wait
        cd $REDUX_BASE_DIR
    done
}

function ratir_copy2web() {
    # copy info files and jpeg's to the web server
    cd $REDUX_BASE_DIR
    dirs="$TODAY bias_bank flat_bank dark_bank"
    $test rsync -avzL --exclude latest_version --exclude 'v*C*' --exclude '*.fits' --exclude '20*.fits.gz' --exclude '20*fits.bz2' --exclude '20*fits.txt' --exclude 'redux_*' --exclude 'start_redux_*' $dirs ${web_user}@${web_pc}:${web_data_dir}/
}

function ratir_do_bias() {
    # set of commands to make bias frames
    source_list=$BIAS_DIRS
    echo "Bias frames: $source_list"
    [ "$source_list" ] && {
        ratir_copy_files
        ratir_create_manifest
        ratir_bias
    }
}

function ratir_do_dark() {
    # set of commands to make dark frames
    source_list=$DARK_DIRS
    echo "Dark frames: $source_list"
    [ "$source_list" ] && {
        ratir_copy_files
        ratir_create_manifest
        ratir_dark
    }
}

function ratir_do_flat() {
    # set of commands to make flat frames
    source_list=$FLAT_DIRS
    echo "Flat frames: $source_list"
    [ "$source_list" ] && {
        ratir_copy_files
        merge_flat_dirs
        source_list=$FLAT_DIRS
        ratir_create_manifest
        ratir_flat
    }
}

function ratir_do_standards() {
    # set of commands to do all standard stars
    source_list=$OPT_STANDARD_DIRS
    echo "Opt Standards: $source_list"
    [ "$source_list" ] && {
        ratir_copy_files
        ratir_create_manifest
        ratir_redux
    }
    source_list=$NIR_STANDARD_DIRS
    echo "NIR Standards: $source_list"
    [ "$source_list" ] && {
        ratir_copy_files
        ratir_create_manifest
        ratir_redux
    }
}

function ratir_full_redux() {
    # do everything
    echo "Doing full redux for TODAY=$TODAY"
    ratir_setdirs
    ratir_do_bias
    ratir_do_dark
    ratir_do_flat
    ratir_do_standards
    ratir_copy2web
}

function ratir_clean_fits() {
    # delete raw fits images
    cd $REDUX_BASE_DIR
    for file in `find . -name '20*fits*'`; do
        $test rm $file
    done
}
