#!/bin/bash

[ -f /home/ddoti/MANUAL_MODE ] && { echo "file /home/ddoti/MANUAL_MODE present, aborting" ; exit 1 ; }

source /usr/local/var/ddoti/bin/redux_funcs_ddoti.sh
test=

# just in case
rm $REDUX_LOCKFILE /home/ddoti/monitoring.txt 2>/dev/null

date -u

ddoti_full_redux
