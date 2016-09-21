#!/bin/bash
#
# Run a binary for a given polynomial order and grid size.
# This is a tiny wrapper around the exahype binary handling the
#  (1) parameters and spec files
#  (2) simulation directory setup
#  (3) actually running
# It is controlled by either env parameters and/or command line
# arguments.
#
# It can run in parallel with lot's of different parameters on the
# same machine.
#
# -- SK, 2016


# PROGRAM PARAMETERS
# Change in code or give as environment variables

# directory where all simulations are stored in
export SIMBASE=${SIMBASE:=simulations/}

# path to the exahype executable to use
export EXABINARY=${EXABINARY:=../../ApplicationExamples/EulerFlow/ExaHyPE-Euler}

# path of the spec file to use
export EXASPECFILE=${EXASPECFILE:=../../ApplicationExamples/EulerFlow.exahype}

function errormsg {
	echo -e "Usage: $0 [-p <num>] [-m <num>]"
	echo -e "   Runs an ExahyPE binary with arguments in the"
	echo -e "   current directory $PWD"
	echo -e "   That is, will do the setup and run."
	exit -1		
}

# bash command line parsing done as in http://stackoverflow.com/a/14203146
while [[ $# -gt 1 ]]; do
	key="$1"
	case $key in
		-p|--polyorder)
		export EXAPORDER="$2"
		shift # past argument
		;;
		-m|--meshsize)
		export EXAMESHSIZE="$2"
		shift # past argument
		;;
		# --default) # example binary string
		#DEFAULT=YES
		#;;
		*)
		# unknown option
		errormsg
	    ;;
	esac
	shift # past argument or value
done

# check if all neccessary arguments are nonempty
[[ ! ${EXAPORDER} ]] && errormsg
[[ ! ${EXAMESHSIZE} ]] && errormsg

# exit script in case of error starting from here
set -e

# compose and setup simulation parameter directory
mkdir -p "$SIMBASE"
SIMDIR="$SIMBASE/p${EXAPORDER}-meshsize${EXAMESHSIZE}/"
mkdir -p "$SIMDIR"

# log everything starting from here
exec > >(awk '{ print strftime("%Y-%m-%d %H:%M:%S"), $0; fflush(); }' |  tee "$SIMDIR/run-$(hostname).log") 2>&1

# convert possibly relative paths to absolute ones
EXABINARY=$(readlink -f "$EXABINARY")
EXASPECFILE=$(readlink -f "$EXASPECFILE")
# and get the pure filenames
BASE_EXABINARY=$(basename "$EXABINARY")
BASE_EXASPECFILE=$(basename "$EXASPECFILE")

# populate simulation directory with absolute symlinks
ln -s "$EXABINARY" "$SIMDIR/"
cp "$EXASPECFILE" "$SIMDIR/"

cd "$SIMDIR"

# setup stuff needed for running the executable
mkdir -p output

# todo: change spec file content
sed -i "s/^\(.*maximum-mesh-size = \).*$/\1${EXAMESHSIZE}/" $BASE_EXASPECFILE
# fun fact: changing the polynomial order requires recompiling. haha.

echo "This is $0 on $(hostname) at $(date)"
echo "Running with the following parameters:"
env | grep -iE 'exa|sim'
echo

# run it
time ./$BASE_EXABINARY $BASE_EXASPECFILE && echo "Finished successfully" || 
	{ echo "ExaHyPE binary failed!"; exit -2; }



