#!/bin/bash

[ -f /home/ddoti/MANUAL_MODE ] && { echo "file /home/ddoti/MANUAL_MODE present, aborting" ; exit 1 ; }

[ -f /home/ddoti/monitoring.txt ] && { echo "GRB monitor already running, exiting" ; exit 2 ; }

source /usr/local/var/ddoti/bin/redux_funcs_ddoti.sh
test=

[ -f "$REDUX_LOCKFILE" ] && { echo "lock file $REDUX_LOCKFILE present, exiting" ; exit 2 ; }

date -u

cd $REDUX_BASE_DIR
ddoti_setdirs

nfiles_min=2

nfiles=0
function getnfiles() {
    local nfiles0=`ls 20*C0o.fits* 2>/dev/null | wc -l`
    local nfiles1=`ls 20*C1o.fits* 2>/dev/null | wc -l`
    local nfiles2=`ls 20*C2o.fits* 2>/dev/null | wc -l`
    local nfiles3=`ls 20*C3o.fits* 2>/dev/null | wc -l`
    local nfiles4=`ls 20*C4o.fits* 2>/dev/null | wc -l`
    local nfiles5=`ls 20*C5o.fits* 2>/dev/null | wc -l`
    local ar=($nfiles0 $nfiles1 $nfiles2 $nfiles3 $nfiles4 $nfiles5)
    local mn=${ar[0]}
    for n in "${ar[@]}" ; do
        ((n < mn && n>0)) && mn=$n
    done
    nfiles=$mn
}

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

source_list=$GRB_DIRS
ready_source_list=

if [ "$source_list" ] ; then
   echo "GRB project initiated for $TODAY!"

   $test ddoti_copy_files

   # find directories ready for processing
   for dir in $source_list ; do
       cd $dir
       getnfiles
       echo "Found >=$nfiles files each of the cameras."
       if [ "$nfiles" -lt "$nfiles_min" ]; then
           echo "Not enough data files yet"
           cd $REDUX_BASE_DIR
           continue
       fi
       [ -s "nfiles_last.txt" ] || echo 0 > nfiles_last.txt
       nfiles_last=`cat nfiles_last.txt`
       if [ "$nfiles" -le "$nfiles_last" ]; then
           echo "No new data files"
           cd $REDUX_BASE_DIR
           continue
       fi
       [ "$test" = "echo" ] || echo $nfiles > nfiles_last.txt
       ready_source_list="$ready_source_list $dir"
       cd $REDUX_BASE_DIR
   done

   touch /home/ddoti/monitoring.txt

   # now process the ready directories
   source_list=$ready_source_list
   [ "$source_list" ] && $test ddoti_do_redux

   rm /home/ddoti/monitoring.txt

fi
