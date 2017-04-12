#!/bin/bash
#
# This is a runscript allowing the templated execution of an ExaHyPE binary.
# Use in context with setupTemplatedSimulation.sh and sim.sh/exa.sh.
#
# -- SvenK for ExaHyPE, Apr 2017

set -e

ME="$(basename "$0")"
err () { >&2 echo $ME: $@; } # print error/info message with name of script
fail () { err $@; exit -1; } # exit after errormessage with name of script
verbose() { echo "$@"; $@; } # execute a command with printing it before

# source the simulation directory parameter file
cd "$(dirname "$0")"
source "runsim.env"

# How much TBB cores to use for shared memory parallelization?
# we assume a default value here if not given
export ExaTbbCores=${ExaTbbCores:=1}

# These two are the most important variables:
[[ ! $BASE_ExaBinary ]] && fail "Need ExaBinary as path to the executable"
[[ ! $BASE_ExaSpecfile ]] && fail "Need ExaSpecfile as path to specification template to use"

# Logging all output as soon as simulation was setup:
#export LOGFILE=${LOGFILE:="$SIMDIR/run-$(hostname).log"}
export LOGFILE_BINARY="run-${BASE_ExaBinary}.log"
export LOGFILE_ALL="run-${ME%.*}.log"

timestampize() { stdbuf -o0 awk '{ print strftime("%Y-%m-%d %H:%M:%S"), $0; fflush(); }'; }
#exec > >(timestampize | tee "$LOGFILE_ALL") 2>&1
tolog() { $@ >> $LOGFILE_ALL 2>&1; } # execute a command and put all output to the extensive log
duplog() { 2>&1 $@ | tee -a $LOGFILE_ALL; } # print to stdout and duplicate to stdout.
info () { tolog echo $@; } # print a message only to the extensive log

# setup stuff needed for running the executable
mkdir -p output vtk-output

# Templated specfile
function exaspecrepl { sed -i "$1" $BASE_ExaSpecfile; }
exaspecrepl '/@comment@/d'
perl -p -i -e 's/@([a-zA-Z0-9]+)@/$ENV{$1}/eig' $BASE_ExaSpecfile

runenv="run.env"
env > $runenv

# allow core dumps
ulimit -c unlimited
# skip any ExaHyPE testing at startup
export EXAHYPE_SKIP_TESTS="True"

info "This is $0 at $(hostname) at $(date)"
info "Starting $BASE_ExaBinary with these compile time specifications:"
info
tolog verbose ./$BASE_ExaBinary --version
info
info "Starting with this specification file:"
info
tolog cat $BASE_ExaSpecfile
info
info "Starting up:"
info
verbose time ./$BASE_ExaBinary $BASE_ExaSpecfile | tee "$LOGFILE_BINARY" | tee -a "$LOGFILE_ALL" \
&& info "Finished ExaHyPE successfully" \
|| {
	err "ExaHyPE binary failed!"
	
	echo "1) You'll find a core dump at $(cat /proc/sys/kernel/core_pattern)"
	echo "2) A verbose log at $LOGFILE_ALL ($(wc -l $LOGFILE_ALL | awk '{print $1'}) lines)"
	echo "3) The short log at $LOGFILE_BINARY ($(wc -l $LOGFILE_ALL | awk '{print $1'}) lines)";
}
