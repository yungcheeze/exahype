#!/bin/bash
#
# Outsourced part of the exa.sh: Updating the repository.
#

SCRIPT="$(readlink -f $0)" # absolute path to update-installation.sh
ME="$(basename "$SCRIPT")" # just my name (update-installation.sh)
SCRIPTDIR="$(dirname $SCRIPT)"
EXA="$SCRIPTDIR/exa.sh"
GITROOT="$(cd $SCRIPTDIR && git rev-parse --show-toplevel)" # absolute path to ExaHyPE repository working copy

verbose() { info $@; $@; } # only used for debugging
info () { echo -e $ME: $@; } # print error/info message with name of script
fail () { info $@; exit -1; } # exit after errormessage with name of script
abort () { echo -e $@; exit -1; } # fail without name of script
finish () { echo $@; exit 0; } # finish with message happily
subreq() { $SCRIPT $@; } # subrequest: Query another command for output
cdroot() { cd "$GITROOT"; } # the crucial change to the repository root directory
getappname() { APPNAME="$PAR"; [ -z "$APPNAME" ] && abort "Usage: $0 $CMD <AppName>"; } # set $APPNAME or die
getapppath() { APPPATH="$(subreq find-appdir "$APPNAME")" || abort "Failure: $APPPATH"; } # set APPPATH or die
cdapp() { cdroot; getappname; getapppath; cd $APPPATH/$APPNAME || abort "Could not go to app"; } # change to application directory

set -e # <- this is important during the script execution!

# no arguments provided
[ $# -eq 0 ] && subreq help

for cmd in "$@"; do
	case $cmd in
		"exahype"|"bootstrap"|"all") # Does git pull, update-peano and update-toolkit
			cdroot; info "Performing standard exahype update"
			git pull
			subreq peano toolkit
			# We don't do libxsmm as it's not needed for most applications
			info "Finished updating ExahyPE (+peano, +toolkit) in $GITROOT"
			;;
		"peano") # Updates the Peano subversion repository
			cdroot; info "Updating Peano"
			cd Peano
			exec ./checkout-update-peano.sh
			info "Finished updating Peano"
			;;
		"toolkit") # Compiles the toolkit with ant and javac
			cdroot; info "Creating Toolkit"
			cd Toolkit
			exec ./build.sh
			;;
		"libxsmm") # Checkout or recompile libxsmm
			cdroot; info "Cloning or updating Libxsmm"
			[[ -e Libxsmm ]] || git clone https://github.com/hfp/libxsmm.git Libxsmm
			# libxsmm is a large repository, instead of cloning, downloading
			# https://github.com/hfp/libxsmm/archive/master.zip
			# would be an option
			cd Libxsmm
			git pull
			make clean
			exec make generator -j4
			;;
		"full-eulerflow") # Update exahype, clean, recompile and run the EulerFlow application
			cdroot; info "Making a full update and cleaning of the standard EulerFlow application"
			# this will be no more neccessary once we have out-of-tree builds.
			StandardApp="EulerFlow"
			subreq exahype
			cd $($EXA find-appdir $StandardApp)
			make clean || true # don't break if never run before
			$EXA compile-run $StandardApp
			;;
		""|"help"|"--help")
			info "Allows updating parts of this repository"
			info "Valid commands are: exahype/bootstrap/all, peano, toolkit, libxsmm"
			info "and all of it's combinations, whitespace seperated"
			;;
		*)
			fail "Could not understand command '$cmd'"
			;;
	esac
done
