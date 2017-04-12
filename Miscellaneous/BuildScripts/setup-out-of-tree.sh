#!/bin/bash
#
# Out-of-tree compilation and Code generation for ExaHyPE
# =======================================================
#
# In ExaHyPE, doing an "out-of-tree" build, ie. having the C files and the object
# files seperated (think of the picture as in http://stackoverflow.com/q/39015453)
# is not as trivial as keeping object files in another directory when you are
# interested in parallel builds with different compile-time constants.
#
# This is due to the ExaHyPE (java) code generator which is heavily used to generate
# different C code for different compile time parameters, for instance the ADERDG
# polynomial order or optimized kernels (python code generator).
#
# That means, while a classical out-of-tree approach would work nicely for classical
# compile time constants (for instance, the MODE={Release,Debug,Asserts} or the
# SHAREDMEM or DISTRIBUTEDMEM settings) it fails as soon as we have to invoke code
# generators.
#
# Following this logic, the next logical step an out-of-tree code generation.
# However, since in the current ExaHyPE build system everything is drilled together
# like hell, this also cannot be jailed into one directory (in principle, however,
# this is easily possible). Our current approach therefore is to mirror (ie. copy)
# all the crucial code directories to an out-of-tree location and build there.
#
# This shell script is written exactly to prepare this: It sets up a copy of the
# most relevant files for compiling in a well-defined manner. It also generates
# scripts in the out-of-tree directory to compile the specific application.
#
# Performance
# ===========
#
# In order to make this not painful as hell, I recommend to out the build directory
# on a fast disc, either a parallel HPC filesystem (lustre, beegefs, etc) or even
# better, a ramdisk.
#
# (c) 2017 ExaHyPE, by SvenK.
#

# Invocation: Just like
#   java -jar Toolkit/dist/ExaHyPE.jar ApplicationExamples/MyCode.exahype
# call this as
#   go/to/compile-out-of-tree.sh ApplicationExamples/MyCode.exahype BuildName
# where BuildName is how you want to call this build.


buildscripts="$(dirname "$0")"
export oot_compile="$buildscripts/compile.sh"

set -e

log() { >&2 echo "$@"; }
warn() { log $@; }
die() { warn "$@"; exit -1; }
verbose() { warn "$@"; >&2 $@; }

[ "$#" -ne 2 ] && die "Usage: $0 <PathToSpecfile.exahype> <BuildName>"
export oot_specfile="$1"
[[ -e "$oot_specfile" ]] || die "Specfile '$oot_specfile' does not exist." 
export oot_base_specfile="$(basename $oot_specfile)" # as a service

# our assumption about the toolkit location
export oot_toolkit="${oot_toolkit=./Toolkit/dist/ExaHyPE.jar}"

# the absolute path to the ExaHyPE base, for the oot scripts
export oot_abs_exahype="$(readlink -f .)"

# other directories (or only subdirectories) we need to copy
export oot_dependencies="CodeGenerator Peano ExaHyPE"
optional_dependencies="Libxsmm/lib ExternalLibraries"

# an assumption: This is called from the base codedir of ExaHyPE,
# ie. the root directory where there are all the dependencies.

for dep in $oot_dependencies; do
	[[ -e $dep ]] || die "Cannot find dependency $dep in pwd $PWD. This maybe because "\
		"you are not in the root directory of ExaHyPE. Please call this script from there."
done

# add optional dependencies if found. This is useful if you don't need Libxsmm.
# however, if you need Libxsmm you have to compile that well before.
for dep in $optional_dependencies; do
	[[ -e $dep ]] && export oot_dependencies="$oot_dependencies $dep"
done


# lightweight parsing the specfile:
# the application directory
export oot_appdir="$(awk '/output-directory/{print substr($0, match($0, /=/)+1); }' $oot_specfile | tr -d '[:space:]')"
# the future name of the binary, once compiled.
export oot_projectname=$(grep '^exahype-project' $oot_specfile | awk '{ print $2; }')
export oot_makesystem_prefix="ExaHyPE-"
export oot_binary="${oot_makesystem_prefix}$oot_projectname"

[[ -d "$oot_appdir" ]] || die "Application directory '$oot_appdir' as given in specfile $oot_specfile does "\
	"not yet exist. This maybe because you're not in the correct root directory (try cd ..)"\
	" or because it's a nonexisting new application. Don't start new projects off-tree."

# Introduce two directories:
#  * outdir:   The managament/master directory where scripts, logs, the final
#              executable goes to
#  * builddir: The directory where the mirror of the gitroot (Peano, ExaHyPE, ...)
#              goes to, ie. all *.c and *.o files.
#              If oot_oot is "yes" (default), this will really go out of tree. Target location
#              can be given with oot_codedirbase which is default /dev/shm, ie. a
#              temporary host-local ramdisk.
#              If oot_oot is "no", it will be just a subdirectory of `outdir`.

export oot_outdirbase="${oot_outdirbase:=./Builds}"
export oot_codedirbase="${oot_codedirbase:=/dev/shm}"
export oot_oot="${oot_oot="Yes"}"

export oot_buildname="$2"
export oot_buildprefix="build-"
export oot_outdir="$oot_outdirbase/$oot_buildprefix$oot_buildname"
mkdir -p $oot_outdir || die "Cannot create $oot_outdir"
log "Setting up build $oot_buildname at $oot_outdir"

export oot_codedir="oot-code" # the subdirectory in $oot_outdir where the files are doubled to
insthash=$(echo $buildscripts | md5sum | fold -w10 | head -n1) # a deterministic hash for this ExaHyPE installation
export qualifiedbuildname="exabuild-$(hostname)-${insthash}-$oot_projectname-$oot_buildname"

if [[ ${oot_oot} == "Yes" ]]; then
	# this is really out of tree
	export oot_builddir="$oot_codedirbase/$qualifiedbuildname"
	mkdir -p $oot_builddir
	[[ -e $oot_outdir/$oot_codedir ]] || ln -s $oot_builddir $oot_outdir/$oot_codedir
	log "Setting up Out-of-Tree copy at $oot_builddir"
else
	# no codedir set, use default
	log "Setting up local (in-tree) copy at $oot_builddir"	
	export oot_builddir="$oot_outdir/$oot_codedir"
	mkdir -p $oot_builddir
fi

##### Create scripts in the OutOfTree directory
#
# All these scripts have access to the oot.env variables, ie.
# all the variables we have defined in this script.

cp $buildscripts/out-of-tree-scripts/*.sh $oot_outdir
chmod 755 $oot_outdir/*.sh

# replace @variables@ in autogenerated scripts with ENV
export AUTOGENINFO="on $(hostname) at $(date) in $(pwd)"

perl -p -i -e 's/@([a-zA-Z0-9]+)@/$ENV{$1}/eg' $oot_outdir/*.sh;
chmod 755 $oot_outdir/*.sh

# used by make.sh
export oot_binarypath="$oot_codedir/$oot_appdir/$oot_binary"

#log "Out of Tree configuration:"

# env will produce files like
# A=b c d
# ie. due to whitespaces, cannot source the file.u
quoteenv() { perl -e 'foreach $k (sort(keys(%ENV))) { print "$k=\"$ENV{$k}\"\n"; }'; }

quoteenv | grep -iE '^oot_' | tee $oot_outdir/oot.env

# now, could call the sync in order to fill the build
# with code:
# $oot_outdir/sync.sh
 
