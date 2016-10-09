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

if [ -z "$ORDER" ]; then
	echo -e "Usage: $0 <polyNomialOrder>"
	exit -1
fi

SCRIPTDIR="$(readlink -f $(dirname "$0"))" # absolute path to script directory
COMPILE="$MYCODE/compile.sh" # absolute path to compile.sh

##################

APPNAME="${PWD##*/}" # eg. "eulerflow2d"
SPECFILE="$APPNAME.exahype" # eg. "eulerflow2d.exahype"
ABSAPPDIR="$(dirname "$PWD")" # absolute path to "Applications"
APPDIRNAME="${ABSAPPDIR##*/}" # eg. "Applications" or "ApplicationExamples"
ABSCODEDIR="$(dirname "$ABSAPPDIR")" # absolute path to "Code"

echo "$0 running with"
echo " APPNAME = $APPNAME"
echo " SPECFILE = $SPECFILE"
echo " SCRIPTDIR = $SCRIPTDIR"
echo " ABSAPPDIR = $ABSAPPDIR"
echo " APPDIRNAME = $APPDIRNAME"
echo " POLYORDERS = $POLYORDERS"
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
BINARYNAME="ExaHyPE-$PROJECTNAME"
FINALBINARYNAME="$BINARYNAME-p$ORDER"
echo "Project name is $PROJECTNAME, so will expect the binary at $BINARYNAME and move to $FINALBINARYNAME"

# and we also have to generate partially our code ourself as the build system
# doesn't do the full job.
echo "Fixing GeneratedConstants.h"
sed -i "s/MY_POLYNOMIAL_DEGREE.*$/MY_POLYNOMIAL_DEGREE = ${ORDER};/" GeneratedConstants.h

export CLEAN="Lightweight"
$COMPILE

# Move the binary to some safe place
echo -e "Copy Binary to safe place"
mv "$BINARYNAME" "$FINALBINARYNAME"

echo -e "Finished at $FINALBINARYNAME"