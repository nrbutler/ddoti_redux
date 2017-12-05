#!/bin/bash

export REDUX_BASE_DIR=/usr/local/var/ddoti
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

# all functions look at a day's worth of data defined by TODAY:
export TODAY=`date -u +20%y%m%d`

# location of raw data and summary (web) files
export raw_data_archive=/nas/archive-ddoti/raw/${TODAY}/executor/images
export web_pc=tcs-a
export web_user=reducer
export web_data_dir=/usr/local/var/www/main/ddoti

nfile_min=2
nflat_min=6
nbias_min=6

export NBATCH=`grep processor /proc/cpuinfo | wc -l`

cd $REDUX_BASE_DIR

function ddoti_setdirs() {
    # default directories, can be changed manually
    $test rsync -a -f"+ */" -f"- *" --chmod="g=rwx" ${raw_data_archive}/ $REDUX_BASE_DIR/${TODAY} 2>/dev/null
    export BIAS_DIRS=`ls -d ${TODAY}/20*-0002/*/* 2>/dev/null`
    export DARK_DIRS=`ls -d ${TODAY}/20*-0003/*/* 2>/dev/null`
    export FLAT_DIRS=`ls -d ${TODAY}/20*-0001/*/* 2>/dev/null`
    export STANDARD_DIRS=`ls -d ${TODAY}/20*-0001/*/* 2>/dev/null`
    export GRB_DIRS=`ls -d ${TODAY}/2017A-1000/*/* 2>/dev/null`
    env | grep _DIRS
}

function ddoti_copy_files() {
    # populate a data directory with raw images from the data server
    for dir in $source_list; do
        dir0=`echo $dir | sed -e "s/${TODAY}\///g"`
        $test ln -s ${raw_data_archive}/$dir0/20*fits.* $dir/ 2>/dev/null
    done
}

function ddoti_create_manifest() {
    # group the raw fits files for reduction
    for dir in $source_list; do
        n=`ls $dir | grep fits | wc -l`
        [ "$n" -eq 0 ] && continue
        $test cd $dir
        $test ls 20*.fits* | awk '{cam=substr($1,16,2); print $1>cam"_list.txt"}'
        $test cd $REDUX_BASE_DIR
    done
}

function ddoti_do_redux() {
    # science data reduction (after bias and flat construction)
    echo "Reducing frames: $source_list"
    ddoti_copy_files
    ddoti_create_manifest
    if [ -f "$REDUX_LOCKFILE" ]; then
        echo "lockfile $REDUX_LOCKFILE present, aborting..."
    else
        $test touch $REDUX_LOCKFILE
        for dir in $source_list; do
            n=`ls $dir | grep fits | wc -l`
            [ "$n" -lt "$nfile_min" ] && continue
            TODAY=`echo $dir | awk -F/ '{print $1}'`
            cd $dir
            for list in `ls C?_list.txt 2>/dev/null`; do
                $test redux_ddoti $list &
            done
            wait
        done
        $test rm $REDUX_LOCKFILE
    fi
}

function ddoti_do_bias() {
    # create bias frames and store to bias bank, work in parallel
    source_list=$BIAS_DIRS
    echo "Bias frames: $source_list"
    ddoti_copy_files
    ddoti_create_manifest
    for dir in $source_list; do
        n=`ls $dir | grep fits | wc -l`
        [ "$n" -lt "$nbias_min" ] && continue
        cd $dir
        for list in `ls C?_list.txt 2>/dev/null`; do
            echo bias_ddoti $list
            $test bias_ddoti $list &
        done
        wait
        cd $REDUX_BASE_DIR
    done
}

function ddoti_do_flat() {
    # create flats frames and store to flat bank, work in parallel
    source_list=$FLAT_DIRS
    echo "Flat frames: $source_list"
    ddoti_copy_files
    ddoti_create_manifest
    for dir in $source_list; do
        n=`ls $dir | grep fits | wc -l`
        [ "$n" -lt "$nflat_min" ] && continue
        cd $dir
        for list in `ls C?_list.txt 2>/dev/null`; do
            echo flat_ddoti $list
            $test flat_ddoti $list &
        done
        wait
        cd $REDUX_BASE_DIR
    done
}

function ddoti_do_standards() {
    # set of commands to do all standard stars
    source_list=$STANDARD_DIRS
    [ "$source_list" ] && ddoti_do_redux
}

function ddoti_full_redux() {
    # do everything
    echo "Doing full redux for TODAY=$TODAY"
    ddoti_do_bias
    ddoti_do_flat
    ddoti_do_standards
}
