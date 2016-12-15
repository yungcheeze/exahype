#!/bin/bash
#
# Creates a number of softlinks to move out heavy data from
# the code directory.
#

# adopt this 
scratchbase="$HOME/numrel/exahype/simulations/ConvergenceTests"

# the directory name in each directory
simulations="simulations"

# do this if you want to clean up directories and softlinks:
#
# rm -f ../*/simulations

verbose () { echo -e $@; $@; }

mkdir -p $scratchbase
verbose ln -s $(readlink -f ..) $scratchbase/codelink

for test in $(dirname $(ls ../*/*.exahype)); do
	btest=$(basename $test)
	echo "Setting up shadow simulation softlink for $btest";
	cd $test
	if [[ -e $simulations ]]; then
		echo "$btest/$simulations exists!"
	else
		verbose mkdir $scratchbase/$btest
		verbose ln -s $scratchbase/$btest $simulations	
		echo "done"
	fi
done
