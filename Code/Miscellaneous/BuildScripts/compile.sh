#!/bin/bash
#
# A compiler script to make the compilation of ExaHyPE applications
# more decent.
#
# (c) 2016 ExaHyPE, Sven K

#cd $(dirname "$0")

# path names for our script
DEFAULT_APPNAME="${PWD##*/}"
DEFAULT_SPECFILE="$DEFAULT_APPNAME.exahype" # eg. "eulerflow2d.exahype"
DEFAULT_ABSAPPDIR="$(dirname "$PWD")" # absolute path to "Applications"
DEFAULT_APPDIRNAME="${DEFAULT_ABSAPPDIR##*/}" # eg. "Applications" or "ApplicationExamples"
DEFAULT_ABSCODEDIR="$(dirname "$DEFAULT_ABSAPPDIR")" # absolute path to "Code"

# options for the Make systems
DEFAULT_COMPILER="GNU"
DEFAULT_SHAREDMEM="None"
DEFAULT_MODE="Asserts"

# options controlling how this script works
DEFAULT_CLEAN="None"
DEFAULT_SKIP_TOOLKIT="No"

# all default variables can be overwritten by specifying them as
# environment variables

APPNAME=${APPNAME:=$DEFAULT_APPNAME}
SPECFILE=${SPECFILE:=$DEFAULT_SPECFILE}
ABSAPPDIR=${ABSAPPDIR:=$DEFAULT_ABSAPPDIR}
APPDIRNAME=${APPDIRNAME:=$DEFAULT_APPDIRNAME}
ABSCODEDIR=${ABSCODEDIR:=$DEFAULT_ABSCODEDIR}
CLEAN=${CLEAN:=$DEFAULT_CLEAN}
SKIP_TOOLKIT=${SKIP_TOOLKIT:=$DEFAULT_SKIP_TOOLKIT}


echo -e "$0 running with"
echo -e " APPNAME = $APPNAME"
echo -e " SPECFILE = $SPECFILE"
echo -e " ABSAPPDIR = $ABSAPPDIR"
echo -e " APPDIRNAME = $APPDIRNAME"
echo -e " CLEAN = $CLEAN"
echo -e " SKIP_TOOLKIT = $SKIP_TOOLKIT"
echo -e "at $(date) on $(hostname) as $(whoami)"
echo -e

export COMPILER=${COMPILER:=$DEFAULT_COMPILER}
export SHAREDMEM=${SHAREDMEM:=$DEFAULT_SHAREDMEM}
export MODE=${MODE:=$DEFAULT_MODE}

# you can amend on this
export TBB_INC=/usr/include/tbb
MPI_LDFLAGS="$(mpicc -showme:link)"
export TBB_SHLIB="-L/usr/lib -ltbb $MPI_LDFLAGS"

set -e

cd "$ABSCODEDIR"

# run the toolkit on this application
if [[ $SKIP_TOOLKIT == "Yes" ]]; then
	echo -e "Skipping toolkit invocation as requested";
elif [ ! -e "$APPDIRNAME/$SPECFILE" ]; then
	echo -e "Cannot find specification file at $APPDIRNAME/$SPECFILE in $PWD";
	exit -1
else
	echo -e "Running toolkit"
	java -jar Toolkit/dist/ExaHyPE.jar  --not-interactive $APPDIRNAME/$APPNAME.exahype
fi

cd -

case $CLEAN in
	"None") echo -e "No cleaning before building."
		;;
	"Clean") echo -e "Make clean"
		make clean
		;;
	"Lightweight") echo -e "Lightweight clean"
		rm *.o
		;;
esac

# Workaround the broken makefile system
echo -e "Fixing Makefile after toolkit run"
sed -i "s/SHAREDMEM=.*/SHAREDMEM=$SHAREDMEM/" Makefile
#sed -i 's/DISTRIBUTEDMEM=.*/DISTRIBUTEDMEM=None/' Makefile

make -j $(nproc)

