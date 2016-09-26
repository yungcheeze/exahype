#!/bin/bash
#
# A compiler script to compile different versions of the ExaHyPE binary.
# A some options in the spec file need a recompilation, this script does it
# for different polynomial orders.
#
# This script does not only make everything compile.sh does (as env setup),
# but also changes the spec file BEFORE running the toolkit and changes
# the Makefile and generated code AFTER running the toolkit and before
# runinng make.
#
# The overall idea is that after running this script, there are several
# almost identical copies of the binary with different names in this
# directory.
#
# (c) 2016 ExaHyPE, Sven K

SCRIPTDIR="$(readlink -f $(dirname "$0"))" # absolute path to script directory, ie. EulerFlow
APPNAME="${PWD##*/}" # eg. "eulerflow2d"
SPECFILE="$APPNAME.exahype" # eg. "eulerflow2d.exahype"
ABSAPPDIR="$(dirname "$PWD")" # absolute path to "Applications"
APPDIRNAME="${ABSAPPDIR##*/}" # eg. "Applications" or "ApplicationExamples"
ABSCODEDIR="$(dirname "$ABSAPPDIR")" # absolute path to "Code"

# Note here the orders you want to have, white space seperated
POLYORDERS="2 3 4 5 6 7 8 9"

echo "$0 running with"
echo " APPNAME = $APPNAME"
echo " SPECFILE = $SPECFILE"
echo " SCRIPTDIR = $SCRIPTDIR"
echo " ABSAPPDIR = $ABSAPPDIR"
echo " APPDIRNAME = $APPDIRNAME"
echo " POLYORDERS = $POLYORDERS"
echo "at $(date) on $(hostname) as $(whoami)"

# for parallelization, see also Makefile editing further below.

export COMPILER=GNU
export SHAREDMEM="None"
#export SHAREDMEM="TBB"
export TBB_INC=/usr/include/tbb
MPI_LDFLAGS="$(mpicc -showme:link)"
export TBB_SHLIB="-L/usr/lib -ltbb $MPI_LDFLAGS"

# Debugging mode
export MODE="Release"

set -e 

for ORDER in $POLYORDERS; do
	cd "$ABSAPPDIR"

	echo "Preparing spec file ${SPECFILE}"
	# change the order parameter to the desire one
	# this is a line with a conent like " order = 3"
	sed -i "s/^\([[:space:]]*order[[:space:]]*=[[:space:]]*\).*$/\1${ORDER}/" ${SPECFILE}
	#echo "checking: "; grep order $SPECFILE;	continue

	# guess name of executable which will be built
	PROJECTNAME=$(grep '^exahype-project' ${SPECFILE} | awk '{ print $2; }')
	BINARYNAME="ExaHyPE-$PROJECTNAME"
	FINALBINARYNAME="$BINARYNAME-p$ORDER"
	echo "Project name is $PROJECTNAME, so will expect the binary at $BINARYNAME and move to $FINALBINARYNAME"

	cd "$ABSCODEDIR"
	# run the toolkit on this application
	if [ ! -e "$APPDIRNAME/$SPECFILE" ]; then
		echo -e "Cannot find specification file at $APPDIRNAME/$SPECFILE in $PWD";
	else
		java -jar Toolkit/dist/ExaHyPE.jar  --not-interactive $APPDIRNAME/$APPNAME.exahype
	fi
	cd "$SCRIPTDIR"

	# the lightweight alternative to "make clean":
	rm  -f *.o

	# Disable any parallel stuff in Makefile. Yeah, the Makesystem is broken
	echo "Fixing Makefile for disabling any parallel stuff. Yeah, the Makesystem is broken"
	sed -i 's/SHAREDMEM=.*/SHAREDMEM=None/' Makefile
	sed -i 's/DISTRIBUTEDMEM=.*/DISTRIBUTEDMEM=None/' Makefile

	# and we also have to generate partially our code ourself as the build system
	# doesn't do the full job.
	echo "Fixing GeneratedConstants.h"
	sed -i "s/MY_POLYNOMIAL_DEGREE.*$/MY_POLYNOMIAL_DEGREE = ${ORDER};/" GeneratedConstants.h

	# Compile!
	make -j $(nproc)

	# Move the binary to some safe place
	echo "Copy Binary to safe place"
	mv "$BINARYNAME" "$FINALBINARYNAME"
done

echo "Finished compiling for all polynomial orders"

