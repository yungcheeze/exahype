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

ORDER="$1"

if [[ ! "$ORDER" ]]; then
	echo -e "Usage: $0 <polyNomialOrder>"
	exit -1
fi

SCRIPTDIR="$(readlink -f $(dirname "$0"))" # absolute path to script directory
COMPILE="$SCRIPTDIR/compile.sh" # absolute path to compile.sh

##################

APPNAME="${PWD##*/}" # eg. "eulerflow2d"
SPECFILE="${SPECFILE:=$APPNAME.exahype}" # eg. "eulerflow2d.exahype"
ABSAPPDIR="$(dirname "$PWD")" # absolute path to "Applications"
APPDIRNAME="${ABSAPPDIR##*/}" # eg. "Applications" or "ApplicationExamples"
ABSCODEDIR="$(dirname "$ABSAPPDIR")" # absolute path to "Code"

echo "$0 running with"
echo " APPNAME = $APPNAME"
echo " SPECFILE = $SPECFILE"
echo " SCRIPTDIR = $SCRIPTDIR"
echo " ABSAPPDIR = $ABSAPPDIR"
echo " APPDIRNAME = $APPDIRNAME"
echo " POLYNOMIAL ORDER = $ORDER"
echo "at $(date) on $(hostname) as $(whoami)"

set -e 

cd "$ABSAPPDIR"

echo "Preparing spec file ${SPECFILE}"
# change the order parameter to the desire one
# this is a line with a conent like " order = 3"
sed -i "s/^\([[:space:]]*order[[:space:]]*=[[:space:]]*\).*$/\1${ORDER}/" ${SPECFILE}
#echo "checking: "; grep order $SPECFILE;	continue

# guess name of executable which will be built
PROJECTNAME=$(grep '^exahype-project' ${SPECFILE} | awk '{ print $2; }')
SOLVERNAME=$(grep 'solver ADER-DG' ${SPECFILE} | awk '{ print $3; }')
BINARYNAME="ExaHyPE-$PROJECTNAME"
FINALBINARYNAME="$BINARYNAME-p$ORDER"
echo "Project name is $PROJECTNAME, so will expect the binary at $BINARYNAME and move to $FINALBINARYNAME"

cd "$APPNAME"

# and we also have to generate partially our code ourself as the build system
# doesn't do the full job.
# MyBlaSolver.h contains the polynomial order
echo "Expect the project solver $SOLVERNAME to be at $SOLVERNAME.h, I delete it for recreation"
rm -f "${SOLVERNAME}.h"

# if no preference for $CLEAN has been set, make sure it is at least Lightweight
# which is needed in order to delete header files etc.
#
# Note: Depending on the Kernels used you might really need a proper CLEAN. Lightweight
# cleaning saves you a lot of compiling time as you don't have to recompile everything
# over and over.
export CLEAN="${CLEAN:=Lightweight}"

$COMPILE || { echo "Failure while compiling p=${ORDER}."; exit -1; }

# Move the binary to some safe place
echo -e "Copy Binary to safe place"
mv "$BINARYNAME" "$FINALBINARYNAME"

# the same with the make log to determine the problem in case of errors
mv make.log make-p${ORDER}.log

echo -e "$0 finished with $FINALBINARYNAME"
