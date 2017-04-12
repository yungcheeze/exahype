#!/bin/bash
#
# sim is a small script for simulation managament.
# uses exa.sh from the same directory.
#
## Caveat: This script is merely a placeholder and will probably
##         implemented using Python instead, or not at all.
#
#
# Note: One day we will implement
#       https://bitbucket.org/dradice/batchtools
#       here for ExaHyPE.
#
# (c) 2016 ExaHyPE - by SvenK

SCRIPT="$(readlink -f $0)" # absolute path to sim.sh
ME="$(basename "$SCRIPT")" # just my name (sim.sh)
SCRIPTDIR="$(dirname $SCRIPT)"
GITROOT="$(cd $SCRIPTDIR && git rev-parse --show-toplevel)" # absolute path to ExaHyPE repository working copy
EXACLITOOL="$SCRIPTDIR/../BuildScripts/exa.sh"
SETUPTOOL="$SCRIPTDIR/setupTemplatedSimulation.sh"

# central parameter for this script:
# directory where simulations are stored.
export SIMULATIONS="${SIMULATIONS:=$(readlink -f "$GITROOT/Simulations")}"

CMD="$1" # the actual command
PAR="$2" # some parameter (for passing to bash functions)
set -- "${@:2}" # pop first parameter

info () { >&2 echo $ME: $@; } # print error/info message with name of script
fail () { info $@; exit -1; } # exit after errormessage with name of script
abort () { >&2 echo $@; exit -1; } # fail without name of script
finish () { echo $@; exit 0; } # finish with message happily
getappname() { APPNAME="$PAR"; [ -z "$APPNAME" ] && abort "Usage: $0 $CMD <AppName>"; } # set $APPNAME or die
subreq() { $0 $@; } # subrequest: Query another command for output
exareq() { $EXACLITOOL $@; } # request to exa.sh
cdroot() { cd "$GITROOT"; } # the crucial change to the repository root directory


if [[ ! -e "$SIMULATIONS" ]]; then
	info "Creating simulation directory $SIMULATIONS"
	mkdir $SIMULATIONS || fail "Could not create directory."
fi

case $CMD in
	"setup-from-app") # create a new simulation quickly from an Application and its specfile
		# This command glues together the "exa find" command with it's standard semantics how and
		# where to find binaries and specfiles with the templated runner scripts in this directory.
		cdroot
		
		# obtain simulation name, ie. name of future simulation directory
		SIMNAME="$1"
		[ -z "$SIMNAME" ] && fail "Usage: sim setup <NameOfSimulation> <NameOfApplication>"
		export SIMDIR="$SIMULATIONS/$SIMNAME"
		
		# obtain application name, ie. name of binary and specfile for the exareq find.
		APPNAME="$2"
		[ -z "$APPNAME" ] && fail "Usage: sim setup <NameOfSimulation> <NameOfApplication>"
		
		# check for application
		export ExaBinary=$(exareq find binary "$APPNAME" 2>&1) || fail "Could not find '$APPNAME': ExaBinary"
		export ExaSpecfile=$(exareq find specfile "$APPNAME" 2>&1) || fail "Error: $ExaSpecfile"
		
		exec $SETUPTOOL
		;;
	"setup-from-build") # create a new simulation quickly from a OOT build and its specfile
		# This command glues together the "exa build" OOT command together with the templated runner
		# scripts from this directory.

		# obtain simulation name, ie. name of the future simulation directory
		SIMNAME="$1"
		[ -z "$SIMNAME" ] && fail "Usage: sim setup <NameOfSimulation> <NameOfBuild>"
		export SIMDIR="$SIMULATIONS/$SIMNAME"

		# obtain application name, ie. name of the OOT build
		BUILDNAME="$2"
		[ -z "$BUILDNAME" ] && fail "Usage: sim setup <NameOfSimulation> <NameOfBuild>"
		
		buildloc=$(subreq root)/$(subreq build-find $BUILDNAME) || fail "Could determine build location of '$BUILDNAME'"
		[[ -e $buildloc/oot.env ]] || abort "Cannot find build instance '$BUILDNAME'. Maybe execute 'exa build-compile $BUILDNAME' before?"
		source $buildloc/oot.env
		
		export ExaBinary="$buildloc/$oot_binary"
		export ExaSpecfile="$buildloc/$oot_base_specfile"
		
		exec $SETUPTOOL
		;;
	"run") # run simulation. Just give it your name
	
		# obtain simulation name, ie. name of the future simulation directory
		SIMNAME="$1"
		[ -z "$SIMNAME" ] && fail "Usage: sim run <NameOfSimulation>"
		SIMDIR="$SIMULATIONS/$SIMNAME"
		
		exec $SIMDIR/ExaRun.sh
		;;
	"setup-run") # Setup and run an application
		fail "Not yet implemented"
		;;
	"help") # Show some help message
		echo "$ME: <command> [parameters]"
		echo "A mini bash'ish simulation managament for ExaHyPE"
		echo
		echo "Available commands:"
		echo 
		cat $0 | grep -oE '^\s+"(.+)"\)' | tr -d '")|'
		echo
		echo "With their individual meanings:"
		echo
		# @TODO: Improve display of available formats
		cat $0 | grep -E '\)\s+#' | tr ')' ':' | tr -d '#' | column -c 2 
		echo
		;;
	*)
		fail "Could not understand command '$CMD'"
		;;
esac

