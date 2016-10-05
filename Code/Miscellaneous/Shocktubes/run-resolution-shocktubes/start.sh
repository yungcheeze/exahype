#!/bin/bash

# invoke with level: 1 or 2 or so

[ -z "$1" ] && { echo -e "Usage: $0 <level>"; exit -1; }

level="$1"
conffile="SRHD-l${level}.exahype"

[ -e $conffile ] || { echo -e "$conffile does not exist"; exit -1; }

exec > >(awk '{ print strftime("%Y-%m-%d %H:%M:%S"), $0; fflush(); }' |  tee run-l${level}-$(hostname).log) 2>&1

echo "This is ExaHyPE starter script at $(hostname) at $(date) with level=$level"
echo

mkdir -p output-l${level}
time ./ExaHyPE-SRHD $conffile

echo Done
