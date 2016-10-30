#!/bin/bash
#
# A compiler script to make the compilation of ExaHyPE applications
# more decent.
#
# (c) 2016 ExaHyPE, Sven K

#cd $(dirname "$0")

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

echo -e "$0 running with"
echo -e " APPNAME = $APPNAME"
echo -e " SPECFILE = $SPECFILE"
echo -e " ABSAPPDIR = $ABSAPPDIR"
echo -e " APPDIRNAME = $APPDIRNAME"
echo -e " CLEAN = $CLEAN"

cd "$ABSCODEDIR"

[[ -e "$APPDIRNAME/$SPECFILE" ]] || { echo -e "Cannot find specfile $APPDIRNAME/$SPECFILE in $PWD"; exit -1; }
PROJECTNAME=$(grep '^exahype-project' "$APPDIRNAME/$SPECFILE" | awk '{ print $2; }')

echo -e " PROJECTNAME = $PROJECTNAME"
echo -e " SKIP_TOOLKIT = $SKIP_TOOLKIT"
echo -e "at $(date) on $(hostname) as $(whoami)"
echo -e

export COMPILER=${COMPILER:=$DEFAULT_COMPILER}
export SHAREDMEM=${SHAREDMEM:=$DEFAULT_SHAREDMEM}
export DISTRIBUTEDMEM=${DISTRIBUTEDMEM:=$DEFAULT_DISTRIBUTEDMEM}
export MODE=${MODE:=$DEFAULT_MODE}

# you can amend on this
export TBB_INC=/usr/include/tbb
MPI_LDFLAGS="$(mpicc -showme:link)"
export TBB_SHLIB="-L/usr/lib -ltbb $MPI_LDFLAGS"

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

	java -jar Toolkit/dist/ExaHyPE.jar  --not-interactive $APPDIRNAME/$APPNAME.exahype
fi

cd -

case $CLEAN in
	"None") echo -e "No cleaning before building."
		;;
	"Clean") echo -e "Make clean"
		make clean
		;;
	"Lightweight") echo -e "Lightweight clean"
		# find also object files in subdirectories
		find . -iname '*.o' -exec rm {} \;
		# and also cleanup Fortran modules
		find . -iname '*.mod' -exec rm {} \;
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

make -j $(nproc) 2>&1 | tee make.log || { echo -e "Make failed, see make.log for full log"; exit -1; }

echo -e "Making $APPNAME finished"

