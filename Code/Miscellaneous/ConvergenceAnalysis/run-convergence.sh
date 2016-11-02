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

export MeaningBinary="Path to the ExaHyPE executable to use"
export ExaBinary=${ExaBinary=$(exa root)/$(exa find-binary MHD)}

export MeaningSpecfile="path of the spec file to use"
export ExaSpecfile=${ExaSpecfile=MHD_AlfenWaveConvergence.exahype}

export MeaningTbbCores="How much TBB cores to use for shared memory parallelization"
export ExaTbbCores=${ExaTbbCores:=16}

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
		export MeaningpOrder="Polynomial order for ADER-DG"
		export ExapOrder="$2"
		ExaBinary="$ExaBinary-p$ExapOrder"
		if [ ! -e "$ExaBinary" ]; then
			echo -e "Failure: $ExaBinary does not exist. Please compile for order $ExapOrder"
			exit -2
		fi
		shift # past argument
		;;
		-m|--meshsize)
		export MeaningMeshSize="Minimum meshsize (to pass to specfile), the coarsest grid resolution"
		export ExaMeshSize="$2"
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
[[ ! ${ExapOrder} ]] && errormsg
[[ ! ${ExaMeshSize} ]] && errormsg

# exit script in case of error starting from here
set -e

# compose and setup simulation parameter directory
mkdir -p "$SIMBASE"
SIMDIR="$SIMBASE/p${ExapOrder}-meshsize${ExaMeshSize}/"
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
ExaBinary=$(readlink -f "$ExaBinary")
ExaSpecfile=$(readlink -f "$ExaSpecfile")
# and get the pure filenames
BASE_ExaBinary=$(basename "$ExaBinary")
BASE_ExaSpecfile=$(basename "$ExaSpecfile")

# populate simulation directory with absolute symlinks
ln -s "$ExaBinary" "$SIMDIR/"
cp "$ExaSpecfile" "$SIMDIR/"

cd "$SIMDIR"

# setup stuff needed for running the executable
mkdir -p output

# set initial data to use.
#export EXAHYPE_INITIALDATA="MovingGauss2D"
export exahype_parameter_workaround="yes"
export exahype_initialdata="AlfenWave"
export EXAHYPE_INITIALDATA="$exahype_initialdata"

# parameters for setting up the specfile
export MeaningSPEC_WIDTH="Grid extend in x direction"
export EXASPEC_WIDTH="1.0"
export MeaningSPEC_HEIGHT="Grid extend in y direction"
export EXASPEC_HEIGHT="1.0"
export EXASPEC_ENDTIME="1.0"
export HOST="$(hostname)"
export DATE="$(date)"

# parameters deciding how frequently output is made. As a first criterion,
# 1 output dump with the highest resolution is 250MB.
export MeaningConvOutputRepeat="Time delta how frequently the data convolution output (integrals) is made"
export ExaConvOutputRepeat="0.05"
export MeaningVtkOutputRepeat="Time delta how frequently the full point-wise VTK output is made"
export ExaVtkOutputRepeat="0.05"

# change spec file contents:
function exaspecrepl { sed -i "$1" $BASE_ExaSpecfile; }

exaspecrepl '/@comment@/d'
exaspecrepl "s/@WIDTH@/$EXASPEC_WIDTH/g"
exaspecrepl "s/@ENDTIME@/$EXASPEC_ENDTIME/g"
exaspecrepl "s/@ORDER@/$ExapOrder/g"
exaspecrepl "s/@MESHSIZE@/$ExaMeshSize/g"
exaspecrepl "s/@HOST@/$HOST/g"
exaspecrepl "s/@DATE@/$DATE/g"
exaspecrepl "s/@TBBCORES@/$ExaTbbCores/g"
exaspecrepl "s/@CONVOUTPUTREPEAT@/$ExaConvOutputRepeat/g"
exaspecrepl "s/@VTKOUTPUTREPEAT@/$ExaVtkOutputRepeat/g"

echo "This is $0 on $HOST at $DATE"
echo "Running with the following parameters:"
env | grep -iE 'exa|sim|meaning' | tee parameters.env
echo

# run it
time ./$BASE_ExaBinary $BASE_ExaSpecfile && echo "Finished ExaHyPE successfully" || 
	{ echo "ExaHyPE binary failed!"; exit -2; }

echo "Packing simulation outcome to results.tar.gzip"
tar cvfz results.tar.gzip *.vtk



