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
export EXABINARY=${EXABINARY:=$(exa root)/$(exa find-binary MHD)}

# path of the spec file to use
export EXASPECFILE=${EXASPECFILE:=MHD_AlfenWaveConvergence.exahype}

# how much TBB cores to use for shared memory parallelization
export EXATBBCORES=${EXATBBCORES:=16}

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
		EXABINARY="$EXABINARY-p$EXAPORDER"
		if [ ! -e "$EXABINARY" ]; then
			echo -e "Failure: $EXABINARY does not exist. Please compile for order $EXAPORDER"
			exit -2
		fi
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
if [ -e "$SIMDIR" ]; then
	echo "WIPING existing simulation at $SIMDIR"
	rm -r "$SIMDIR";
fi
mkdir -p "$SIMDIR"

# log everything starting from here
# exec > >(awk '{ print strftime("%Y-%m-%d %H:%M:%S"), $0; fflush(); }' |  tee "$SIMDIR/run-$(hostname).log") 2>&1
#
# or instead, ONLY put everything to the log starting from here, surpressing stdout.

# stdbuf -o0 to avoid 4kB buffering
LOGFILE="$SIMDIR/run-$(hostname).log"
echo "Redirecting further output to $LOGFILE and executing ExaHyPE in background"
exec > >(stdbuf -o0 awk '{ print strftime("%Y-%m-%d %H:%M:%S"), $0; fflush(); }' > $LOGFILE) 2>&1

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

# set initial data to use.
#export EXAHYPE_INITIALDATA="MovingGauss2D"
export EXAHYPE_INITIALDATA="ShuVortex"

# parameters for setting up the specfile
export EXASPEC_WIDTH="1.0"
export EXASPEC_HEIGHT="1.0"
export EXASPEC_ENDTIME="1.0"
export HOST="$(hostname)"
export DATE="$(date)"

# parameters deciding how frequently output is made. As a first criterion,
# 1 output dump with the highest resolution is 250MB.
export EXACONVOUTPUTREPEAT="0.05"
export EXAVTKOUTPUTREPEAT="0.05"

# change spec file contents:
function exaspecrepl { sed -i "$1" $BASE_EXASPECFILE; }

exaspecrepl '/@comment@/d'
exaspecrepl "s/@WIDTH@/$EXASPEC_WIDTH/g"
exaspecrepl "s/@ENDTIME@/$EXASPEC_ENDTIME/g"
exaspecrepl "s/@ORDER@/$EXAPORDER/g"
exaspecrepl "s/@MESHSIZE@/$EXAMESHSIZE/g"
exaspecrepl "s/@HOST@/$HOST/g"
exaspecrepl "s/@DATE@/$DATE/g"
exaspecrepl "s/@TBBCORES@/$EXATBBCORES/g"
exaspecrepl "s/@CONVOUTPUTREPEAT@/$EXACONVOUTPUTREPEAT/g"
exaspecrepl "s/@VTKOUTPUTREPEAT@/$EXAVTKOUTPUTREPEAT/g"

echo "This is $0 on $HOST at $DATE"
echo "Running with the following parameters:"
env | grep -iE 'exa|sim' | tee parameters.env
echo

# run it
time ./$BASE_EXABINARY $BASE_EXASPECFILE && echo "Finished ExaHyPE successfully" || 
	{ echo "ExaHyPE binary failed!"; exit -2; }

echo "Packing simulation outcome to results.tar.gzip"
tar cvfz results.tar.gzip *.vtk



