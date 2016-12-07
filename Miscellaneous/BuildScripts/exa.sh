#!/bin/bash
#
# exa is a versatile CLI entry point for various scripts.
# Install it with a softlink from your /usr/local/bin or ~/bin to
# make use of the full delocalized power.
#
#
# Idee: Ebenfalls als bashrc script nicht nur fuer completion sondern auch slimmeres gefuehl
#       beim Tabben und und "exa cd" oder "exa paths" erlaubt setzen von umgebungsvariablen direkt
#       drinnen.
#       Alternativ: exa shell mit Parametern und veraenderter bashrc? finde ich weniger sinnvoll.
#
# (c) 2016 ExaHyPE - by SvenK

SCRIPT="$(readlink -f $0)" # absolute path to exa.sh
ME="$(basename "$SCRIPT")" # just my name (exa.sh)
SCRIPTDIR="$(dirname $SCRIPT)"
GITROOT="$(cd $SCRIPTDIR && git rev-parse --show-toplevel)" # absolute path to ExaHyPE repository working copy

CMD="$1" # the actual command
PAR="$2" # some parameter (for passing to bash functions)
set -- "${@:2}" # pop first parameter

info () { echo -e $ME: $@; } # print error/info message with name of script
fail () { info $@; exit -1; } # exit after errormessage with name of script
abort () { echo -e $@; exit -1; } # fail without name of script
finish () { echo $@; exit 0; } # finish with message happily
subreq() { $0 $@; } # subrequest: Query another command for output
cdroot() { cd "$GITROOT"; } # the crucial change to the repository root directory
getappname() { APPNAME="$PAR"; [ -z "$APPNAME" ] && abort "Usage: $0 $CMD <AppName>"; } # set $APPNAME or die
getapppath() { APPPATH="$(subreq find-appdir "$APPNAME")" || abort "Failure: $APPPATH"; } # set APPPATH or die
cdapp() { cdroot; getappname; getapppath; cd $APPPATH/$APPNAME || abort "Could not go to app"; } # change to application directory

case $CMD in
	"update-peano") # Updates the Peano subversion repository
		cdroot; info "Updating Peano"
		cd Peano
		exec ./checkout-update-peano.sh
		info "Finished updating Peano"
		;;
	"update-toolkit") # Compiles the toolkit with ant and javac
		cdroot; info "Creating Toolkit"
		cd Toolkit
		exec ./build.sh
		;;
	"update-libxsmm") # Checkout or recompile libxsmm
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
	"list-apps") # Lists all ExaHyPE applications available. Use "find-app" for full path.
		cdroot; info "Listing available Applications:"
		find ApplicationExamples/* -type d -exec basename {} \;
		;;
	"find-app") # Gives the full path from ExaHyPE root to an application
		# this is trivial since we have no more {Applications, ApplicationExamples} directories
		cdroot; getappname
		ls -d ApplicationExamples/$APPNAME || fail "Application '$APPNAME' not found"
		;;
	"find-specfile") # Gives the full path from ExaHyPE root to an application specfile
		# this is trivial since we have no more {Applications, ApplicationExamples} directories
		cdroot; getappname
		ls -f ApplicationExamples/$APPNAME.exahype || fail "Application '$APPNAME' not found";
		exit 0
		;;
	"find-appdir") # Gives directory where app lives inside
		cdroot; getappname
		ls -d Applications/$APPNAME &>/dev/null && finish "Applications/"
		ls -d ApplicationExamples/$APPNAME &>/dev/null && finish "ApplicationExamples/"
		fail "Could not find Application '$APPNAME' somewhere"
		;;
	"find-binary") # Gives path to the executable, even if not present
		cdroot; getappname; getapppath
		SPECFILE="$(subreq find-specfile "$APPNAME")" || abort "Specfile Failure: $SPECFILE"
		PROJECTNAME=$(grep '^exahype-project' ${SPECFILE} | awk '{ print $2; }')
		echo $APPPATH/$APPNAME/ExaHyPE-$PROJECTNAME
		;;
	"toolkit") # Run the toolkit for an application, without compiling
		cdroot; getappname
		SPECFILE="$(subreq find-specfile "$APPNAME")" || abort "Could not find specfile: $SPECFILE"
		info "Running ExaHyPE.jar on $SPECFILE"
		java -jar Toolkit/dist/ExaHyPE.jar --not-interactive $SPECFILE
		;;
	"compile") # Invokes the toolkit and compilation of an application
		cdapp; $SCRIPTDIR/compile.sh
		;;
	"polycompile") # Compile for different polynomial orders (as basis for convergence studies).
		set -- "${@:2}" # pop parameter
		cdapp; $SCRIPTDIR/compile-for-polyorder.sh $@
		;;
	"make") # compiel without invoking the toolkit
		cdapp
		export SKIP_TOOLKIT="Yes"
		export CLEAN="${CLEAN:=Lightweight}" # do no heavy cleaning
		$SCRIPTDIR/compile.sh
		;;
	"git") # passes commands to git
		cdroot; info "ExaHyPE Git Repository at $GITROOT"
		exec git $@
		;;
	"pwd") # give current working directory relative to ExaHyPE root	
		echo "Your directory:   $PWD"
		echo "ExaHyPE root dir: $GITROOT"
		fail "Calculation not yet implemented"
		;;
	"player") # passes commands to the plotting toolkit exaplot. Use "--help" for help.
		exec $GITROOT/Miscellaneous/Postprocessing/exaplayer.py $@
		;;
	"reader") # passes commands to the exahype python conversion toolkit. Use "--help" for help.
		exec $GITROOT/Miscellaneous/Postprocessing/exareader.py $@
		;;
	"run") # quickly start an application inside it's directory. Cleans VTK files before.
		cdapp
		info "Starting $APPNAME with SKIP_TESTS." | tee run.log
		rm -f *.vtk *.log-file
		export EXAHYPE_SKIP_TESTS=TRUE
		ROOT=$(subreq root)
		BINARY=$(subreq find-binary $APPNAME) && [[ -e "$ROOT/$BINARY" ]] || fail "Could not find binary ($BINARY)"
		SPECFILE=$(subreq find-specfile $APPNAME)
		# the directory handling of this tool is really awkward.
		# Would prefer relative directories here.
		reducedbuf="stdbuf -oL -eL" # for quicker output, no 4k buffering
		$reducedbuf $ROOT/$BINARY $ROOT/$SPECFILE 2>&1 | $reducedbuf tee -a run.log
		;;
	"sim") # lightweight simulation managament
		exec $GITROOT/Miscellaneous/BuildScripts/sim.sh $@
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
		cat $0 | grep -E '\)\s+#' | tr ')' ':' | tr -d '#' | column -c 2 
		echo -e
		;;
	"is") # prints out fortunes
		echo "cool"
		;;
	"todo") # prints out all files where TODO notes are inside
		cdroot; info "Probably all current files containing todo messages"
		find . -type f | grep -Ei '\.(cpp|C|h|f90|cc|tex)$' | xargs grep -i todo
		;;
	"root") # prints out the root of the ExaHyPE installation
		echo $GITROOT
		;;
	*)
		fail "Could not understand command '$CMD'"
		;;
esac

