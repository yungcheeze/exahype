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

# Note here the orders you want to have, white space seperated
POLYORDERS="2 3 4 5 6 7 8 9"

echo -e "$0 running with"
echo -e " POLYORDERS = $POLYORDERS"
echo -e "at $(date) on $(hostname) as $(whoami)"

set -e 

for ORDER in $POLYORDERS; do
	./compile-for-polyorder.sh $ORDER
done

echo -e "Finished compiling for all polynomial orders"

