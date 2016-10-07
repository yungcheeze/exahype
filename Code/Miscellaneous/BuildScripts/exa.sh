#!/bin/bash
#
# exa is a versatile CLI entry point for various scripts.
#
# Idee: Ebenfalls als bashrc script nicht nur fuer completion sondern auch slimmeres gefuehl
# beim Tabben und und "exa cd" oder "exa paths" erlaubt setzen von umgebungsvariablen direkt
# drinnen.
# Alternativ: exa shell mit Parametern und veraenderter bashrc? finde ich weniger sinnvoll.
#
#
#
#
#
#

SCRIPT="$(readlink -f $0)" # absolute path to exa.sh
ME="$(basename "$SCRIPT")" # just my name (exa.sh)
SCRIPTDIR="$(dirname $SCRIPT)"
GITROOT="$(cd $SCRIPTDIR && git rev-parse --show-toplevel)" # absolute path to ExaHyPE repository working copy

# this is crucial: exa.sh operates from the repository root directory
cd "$GITROOT"

CMD="$1" # the actual command
PAR="$2" # some parameter (for passing to bash functions)
set -- "${@:2}" # pop first parameter

info () { echo -e $ME: $@; } # print error/info message with name of script
fail () { info $@; exit -1; } # exit after errormessage with name of script
abort () { echo -e $@; exit -1; } # fail without name of script
finish () { echo $@; exit 0; } # finish with message happily
getappname() { APPNAME="$PAR"; [ -z "$APPNAME" ] && abort "Usage: $0 $CMD <AppName>"; } # set $APPNAME or die
subreq() { $0 $@; } # subrequest: Query another command for output

case $CMD in
	"update-peano") # Unpacks Peano from tarball
		info "Updating Peano"
		cd Code/Peano
		tar xvfz peano.tar.gz
		git checkout .gitignore
		info "Finished updating Peano"
		;;
	"list-apps") # Lists all ExaHyPE applications available. Use "find-app" for full path.
		info "Listing available Applications:"
		# we merge directories "Applications" and "ApplicationExamples" in the output
		find Code/Application{,Example}s/* -type d -exec basename {} \;
		;;
	"find-app") # Gives the full path from ExaHyPE root to an application
		getappname
		ls -d Code/Applications/$APPNAME 2>/dev/null; f1=$?
		ls -d Code/ApplicationExamples/$APPNAME 2>/dev/null; f2=$?
		[[ (( $f1 == 0 )) && (( $f2 == 0 )) ]] && fail "Clash of Application Names"
		[[ (( $f1 != 0 )) && (( $f2 != 0 )) ]] && fail "Application '$APPNAME' not found"
		;;
	"find-specfile") # Gives the full path from ExaHyPE root to an application specfile
		getappname
		ls -f Code/Applications/$APPNAME.exahype 2>/dev/null; f1=$?
		ls -f Code/ApplicationExamples/$APPNAME.exahype 2>/dev/null; f2=$?
		[[ (( $f1 == 0 )) && (( $f2 == 0 )) ]] && fail "Clash of Application Names";
		[[ (( $f1 != 0 )) && (( $f2 != 0 )) ]] && fail "Application '$APPNAME' not found";
		exit 0
		;;
	"find-appdir") # Gives directory where app lives inside
		getappname
		ls -d Code/Applications/$APPNAME &>/dev/null && finish "Code/Applications/"
		ls -d Code/ApplicationExamples/$APPNAME &>/dev/null && finish "Code/ApplicationExamples/"
		fail "Could not find Application '$APPNAME' somewhere"
		;;
	"toolkit") # Run the toolkit for an application, without compiling
		getappname
		SPECFILE="$(subreq find-specfile "$APPNAME")" || abort "Could not find specfile: $SPECFILE"
		info "Running ExaHyPE.jar on $SPECFILE"
		cd Code
		java -jar Toolkit/dist/ExaHyPE.jar --not-interactive ../$SPECFILE
		;;
	"compile") # Invokes the toolkit and compilation of an application
		getappname
		APPPATH="$(subreq find-appdir "$APPNAME")" || abort "Failure: $APPPATH"
		cd $APPPATH
		$SCRIPTDIR/compile.sh
		;;
	"git") # passes commands to git
		info "ExaHyPE Git Repository at $GITROOT"
		exec git $@
		;;
	"pwd") # give current working directory relative to ExaHyPE root	
		# todo: Need knowledge of old PWD, thus do not change on beginning.
		fail "Not yet implemented"
		;;
	"player") # passes commands to the plotting toolkit exaplot. Use "--help" for help.
		exec ./Code/Miscellaneous/Postprocessing/exaplayer.py $@
		;;
	"reader") # passes commands to the exahype python conversion toolkit. Use "--help" for help.
		exec ./Code/Miscellaneous/Postprocessing/exareader.py $@
		;;
	""|"help") # prints out help about the exa toolkit
		me=$(basename "$0")
		echo -e "$me: <command> [parameters]"
		echo -e "An ExaHyPE quick command helper."
		echo -e "It is operating with the ExaHyPE installation at $GITROOT"
		echo -e
		echo -e "Available commands:"
		echo -e
		cat $0 | grep -oE '^\s+"(.+)"\)' | tr -d '")|'
		echo -e
		echo -e "With their individual meanings:"
		echo -e
		# @TODO: Improve display of available formats
		cat $0 | grep -E '\)\s+#' | tr ')' ':' | tr -d '#' | columns -c 2 
		echo -e
		;;
	"is") # prints out fortunes
		echo "cool"
		;;
	"todo") # prints out all files where TODO notes are inside
		info "Probably all current files containing todo messages"
		find . -type f | grep -Ei '\.(cpp|C|h|f90|cc|tex)$' | xargs grep -i todo
		;;
	"root") # prints out the root of the ExaHyPE installation
		echo $GITROOT
		;;
	*)
		fail "Could not understand command '$CMD'"
		;;
esac
