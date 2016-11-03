#!/bin/bash
#
# sim is a small script for simulation managament.
# uses exa.sh from the same directory.
#
## Caveat: This script is merely a placeholder and will probably
##         implemented using Python instead, or not at all.
#
# (c) 2016 ExaHyPE - by SvenK

SCRIPT="$(readlink -f $0)" # absolute path to sim.sh
ME="$(basename "$SCRIPT")" # just my name (sim.sh)
SCRIPTDIR="$(dirname $SCRIPT)"
GITROOT="$(cd $SCRIPTDIR && git rev-parse --show-toplevel)" # absolute path to ExaHyPE repository working copy
EXACLITOOL="$SCRIPTDIR/exa.sh"

# central parameter for this script:
# directory where simulations are stored.
SIMULATIONS="$(readlink -f "$GITROOT/../")"
echo $SIMULATIONS

CMD="$1" # the actual command
PAR="$2" # some parameter (for passing to bash functions)
set -- "${@:2}" # pop first parameter

info () { echo -e $ME: $@; } # print error/info message with name of script
fail () { info $@; exit -1; } # exit after errormessage with name of script
abort () { echo -e $@; exit -1; } # fail without name of script
finish () { echo $@; exit 0; } # finish with message happily
getappname() { APPNAME="$PAR"; [ -z "$APPNAME" ] && abort "Usage: $0 $CMD <AppName>"; } # set $APPNAME or die
subreq() { $0 $@; } # subrequest: Query another command for output
exareq() { $EXA $@; } # request to exa.sh

case $CMD in
	"setup") # create a new simulation
		# usage like: exa sim setup <NameOfSimulation> --application=EulerFlow [--specfile=../.../foo.exahype]
		fail "Not yet implemented"
		;;
	"setup-range") # setup a range of simulations
		# usage like exa sim setup-range <NameOfBase> --application=EulerFlow
		# However, ranges are currently done with individual scripts
		fail "Not yet implemented"
		;;
	"run") # run simulation or range of simulations
		fail "Not yet implemented"
		;;
	"")
		fail "Not yet implemented"
		;;
	*)
		fail "Could not understand command '$CMD'"
		;;
esac

