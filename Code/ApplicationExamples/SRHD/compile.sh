#!/bin/bash
#
# A compiler script to make the compilation of ExaHyPE applications
# more decent.

cd $(dirname "$0")
APPNAME="${PWD##*/}" # eg. "eulerflow2d"
SPECFILE="$APPNAME.exahype" # eg. "eulerflow2d.exahype"
ABSAPPDIR="$(dirname "$PWD")" # absolute path to "Applications"
APPDIRNAME="${ABSAPPDIR##*/}" # eg. "Applications" or "ApplicationExamples"
ABSCODEDIR="$(dirname "$ABSAPPDIR")" # absolute path to "Code"

echo "Compile.sh running with"
echo " APPNAME = $APPNAME"
echo " SPECFILE = $SPECFILE"
echo " ABSAPPDIR = $ABSAPPDIR"
echo " APPDIRNAME = $APPDIRNAME"
echo "at $(date) on $(hostname) as $(whoami)"
echo

export COMPILER=GNU
export SHAREDMEM="TBB" # None
export TBB_INC=/usr/include/tbb
MPI_LDFLAGS="$(mpicc -showme:link)"
export TBB_SHLIB="-L/usr/lib -ltbb $MPI_LDFLAGS"

# Debugging mode
export MODE="DEBUG"

set -e

cd "$ABSCODEDIR" # equal to cd ../../

# run the toolkit on this application
if [ ! -e "$APPDIRNAME/$SPECFILE" ]; then
	echo -e "Cannot find specification file at $APPDIRNAME/$SPECFILE in $PWD";
else
	java -jar Toolkit/dist/ExaHyPE.jar  --not-interactive $APPDIRNAME/$APPNAME.exahype
fi

cd -

#make clean
# a lightweight alternative to "make clean" is
# rm *.o

make -j $(nproc)

