#!/bin/bash
#
# This script is a tiny simulation manager which runs several simulations in parallel
# and keeping code and simulation data in seperated directories.
# Use it on a big machine like
# 
#    for p in 2 3 4 5 6 7 8 9; do ./run-shocktubes.sh $p & done;
#
# The script expects one parameter which is the $ORDER.
#
# In order to be able to make use of this script, you need to compile several versions
# of the specific application (ie. SRHD) in parallel. Use the script
# "multicompile-polyorders.sh" to do that.
#
# -- Sven, Oct 2016 for ExaHyPE

ORDER="$1"

if [ -z "$ORDER" ]; then
	echo -e "Usage: $0 <polyNomialOrder>"
	exit -1
fi

set -e # exit on error

# ORDER is something like in
# POLYORDERS="2 3 4 5 6 7 8 9"

ApplicationExamples="../../ApplicationExamples/"
AbsApplicationsPath=$(readlink -f "$ApplicationExamples")
Application="SRHD"
Simulations="/home/koeppel/numrel/exahype/simulations/ADERDG-Shocktubes"

mkdir -vp $Simulations

Simdir="$Application-p$ORDER"
echo "Preparing simulation directory $Simdir"
if [ -e "$Simulations/$Simdir" ]; then
	echo "Wikiping existing simulations directory!"
	rm -rf "$Simulations/$Simdir"
fi
mkdir -p "$Simulations/$Simdir"

BINARYNAME="ExaHyPE-$Application-p$ORDER"
echo "Softllinking executable $BINARYNAME"
cd $Simulations/$Simdir
ln -s $AbsApplicationsPath/$Application/$BINARYNAME

SPECFILE="$Application.exahype"
echo "Preparing spec file ${SPECFILE}"
# change the order parameter to the desire one
# this is a line with a conent like " order = 3"
sed "s/^\([[:space:]]*order[[:space:]]*=[[:space:]]*\).*$/\1${ORDER}/" $AbsApplicationsPath/$SPECFILE > ./$SPECFILE

echo "Setup stuff needed for running the executable"

mkdir -p output

LOGFILE="run-$(hostname).log"
echo "Redirecting further output to $LOGFILE"
exec > >(stdbuf -o0 awk '{ print strftime("%Y-%m-%d %H:%M:%S"), $0; fflush(); }' > $LOGFILE) 2>&1

echo "This is $0 on $(hostname) at $(date)"
echo "Running the polynomial orders Shocktube test"
echo

time ./$BINARYNAME $SPECFILE || { echo "ExaHyPE binary failed!"; exit -2; }

echo "Finished successfully"
