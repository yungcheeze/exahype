#!/bin/bash
#
# A wrapper around typing "make"  inside an ExaHyPE application
# directory, featuring
#
# 1. environment parameter driven behaviour (just like ExaHyPE's Makefile)
# 2. Invocation of the toolkit
# 3. Patching of faulty files where the toolkit is broken
# 4. Cleaning (also partially) before make
# 5. Logging the output of make, also passing to whoopsie in case of failure
# 6. Invoking auto parallelized make.
#
# To use, go to the ExaHyPE app you want to compile, and instead of typing
#   > make -j4 2&>1 | tee make.log
# Just type
#   > ../path/to/buildScripts/compile.sh
# Or just use exa:
#   > exa compile YourApp
# from anywhere.
#
# (c) 2016 ExaHyPE, Sven K

buildscripts="$(dirname "$0")"

verbose() { echo -e $@; $@; }

# path names for our script
DEFAULT_APPNAME="${PWD##*/}"
DEFAULT_SPECFILE="$DEFAULT_APPNAME.exahype" # eg. "eulerflow2d.exahype"
DEFAULT_ABSAPPDIR="$(dirname "$PWD")" # absolute path to "Applications"
DEFAULT_APPDIRNAME="${DEFAULT_ABSAPPDIR##*/}" # eg. "Applications" or "ApplicationExamples"
DEFAULT_ABSCODEDIR="$(dirname "$DEFAULT_ABSAPPDIR")" # absolute path to "Code"

# options for the Make systems
DEFAULT_COMPILER="GNU"
DEFAULT_SHAREDMEM="None"
DEFAULT_DISTRIBUTEDMEM="None"
DEFAULT_MODE="Asserts"

# options controlling how this script works
DEFAULT_CLEAN="None"
DEFAULT_SKIP_TOOLKIT="No"

# all default variables can be overwritten by specifying them as
# environment variables

APPNAME=${APPNAME:=$DEFAULT_APPNAME}
SPECFILE=${SPECFILE:=$DEFAULT_SPECFILE}
ABSAPPDIR=${ABSAPPDIR:=$DEFAULT_ABSAPPDIR}
APPDIRNAME=${APPDIRNAME:=$DEFAULT_APPDIRNAME}
ABSCODEDIR=${ABSCODEDIR:=$DEFAULT_ABSCODEDIR}
CLEAN=${CLEAN:=$DEFAULT_CLEAN}
SKIP_TOOLKIT=${SKIP_TOOLKIT:=$DEFAULT_SKIP_TOOLKIT}

# go to ExaHyPE-Engine root directory (used to be Code/ in former days)
verbose cd "$ABSCODEDIR" || { echo -e "Cannot compile $APPNAME as there is no ABSCODEDIR=$ABSCODEDIR"; exit -1; }

# Logging all further invocations of the toolkit, etc. to make.log
unbuf="stdbuf -i0 -o0 -e0" # turn off buffering in pipe
exec &> >($unbuf tee "$ABSAPPDIR/$APPNAME/make.log")

echo -e "$0 running with"
echo -e " APPNAME = $APPNAME"
echo -e " SPECFILE = $SPECFILE"
echo -e " ABSAPPDIR = $ABSAPPDIR"
echo -e " ABSCODEDIR = $ABSCODEDIR"
echo -e " APPDIRNAME = $APPDIRNAME"
echo -e " CLEAN = $CLEAN"

[[ -e "$APPDIRNAME/$SPECFILE" ]] || { echo -e "Cannot find specfile $APPDIRNAME/$SPECFILE in $PWD"; exit -1; }
PROJECTNAME=$(grep '^exahype-project' "$APPDIRNAME/$SPECFILE" | awk '{ print $2; }')

echo -e " PROJECTNAME = $PROJECTNAME"
echo -e " SKIP_TOOLKIT = $SKIP_TOOLKIT"

export COMPILER=${COMPILER:=$DEFAULT_COMPILER}
export SHAREDMEM=${SHAREDMEM:=$DEFAULT_SHAREDMEM}
export DISTRIBUTEDMEM=${DISTRIBUTEDMEM:=$DEFAULT_DISTRIBUTEDMEM}
export MODE=${MODE:=$DEFAULT_MODE}

# you can amend on this
export TBB_INC=/usr/include/tbb
MPI_LDFLAGS="$(mpicc -showme:link)"
export TBB_SHLIB="-L/usr/lib -ltbb $MPI_LDFLAGS"

#echo -e " COMPILER=$COMPILER"
#echo -e " SHAREDMEM=$SHAREDMEM"
#echo -e " DISTRIBUTEDMEM=$DISTRIBUTEDMEM"
#echo -e " MODE=$MODE"

echo -e "at $(date) on $(hostname) as $(whoami)"
echo -e

set -e

# run the toolkit on this application
if [[ $SKIP_TOOLKIT == "Yes" ]]; then
	echo -e "Skipping toolkit invocation as requested";
else
	echo -e "Running toolkit"

	# todo: 
	#echo -e "Working around defect Makefiles etc"
	#rm $APPDIRNAME/$APPNAME/Makefile
	#could also delete KernelCalls.cpp, $APPNAME_generated.cpp, etc.

	verbose java -jar Toolkit/dist/ExaHyPE.jar  --not-interactive $APPDIRNAME/$APPNAME.exahype
fi

cd -

# plausability check
[[ -e Makefile ]] || { echo -e "Could not find Makefile in $PWD. Probably the toolkit failed."; exit -1; }

case $CLEAN in
	"None") echo -e "No cleaning before building."
		;;
	"Clean") 
		verbose make clean
		;;
	"Lightweight") echo -e "Lightweight clean"
		# find also object files in subdirectories
		verbose find . -iname '.o' -exec rm {} \;
		# and also cleanup Fortran modules
		verbose find . -iname '.mod' -exec rm {} \;
		;;
esac

# Workaround the broken makefile system
echo -e "Fixing Makefile after toolkit run"
sed -i "s/SHAREDMEM=.*/SHAREDMEM=$SHAREDMEM/" Makefile
sed -i "s/DISTRIBUTEDMEM=.*/DISTRIBUTEDMEM=$DISTRIBUTEDMEM/" Makefile

# Workaround for files overwritten by toolkit
for patchfile in "${PROJECTNAME}_generated.cpp"; do # add files seperated by whitespace as in "foo" "bar"
	if [ -e $patchfile ]; then
		if git ls-files $patchfile --error-unmatch &>/dev/null; then
			echo -e "Patching $patchfile with committed version to overwrite toolkit changes"
			git checkout $patchfile
		else
			echo -e "Patching: File $patchfile not under version control"
		fi
	else
		echo -e "Don't patch $patchfile as not there"
	fi
done

# Workaround for broken dependency build system:
# Generate well known fortran modules if present

# This works, but cannot access the EXAHYPE_CFLAGS or similar. Instead, I made
# a small change in the Makesystem so Modules compile first. Probably.

#for fmodule in Parameters.f90 typesDef.f90; do
#	if [ -e $fmodule ]; then
#		echo -e "Precompiling $fmodule as otherwise build fails"
#		FORTFLAGS="-fdefault-real-8 -fdefault-double-8 -ffree-line-length-none"
#		verbose gfortran $FORTFLAGS -c $fmodule
#	fi
#done

set -o pipefail # fail if make fails

verbose make -j $(nproc) || {
	echo -e "Make failed!";
	$buildscripts/whoopsie-paster.sh
	exit -1;
}

echo -e "Making $APPNAME finished successfully"

