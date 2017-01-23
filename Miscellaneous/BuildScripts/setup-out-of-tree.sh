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
compile="$buildscripts/compile.sh"

set -e

log() { >&2 echo "$@"; }
warn() { log $@; }
die() { warn "$@"; exit -1; }
verbose() { warn "$@"; >&2 $@; }

[ "$#" -ne 2 ] && die "Usage: $0 <PathToSpecfile.exahype> <BuildName>"
export specfile="$1"
[[ -e "$specfile" ]] || die "Specfile '$specfile' does not exist." 

# our assumption about the toolkit
export TOOLKIT="${TOOLKIT:=./Toolkit/dist/ExaHyPE.jar}"

# other directories (or only subdirectories) we need to copy
export DEPENDENCIES="CodeGenerator Libxsmm/lib Peano ExaHyPE"

# lightweight parsing the specfile:
# the application directory
export oot_appdir="$(awk '/output-directory/{print substr($0, match($0, /=/)+1); }' $specfile | tr -d '[:space:]')"
# the future name of the binary, once compiled.
export oot_projectname=$(grep '^exahype-project' $specfile | awk '{ print $2; }')
export makesystem_prefix="ExaHyPE-"
export oot_binary="${makesystem_prefix}$oot_projectname"

[[ -d "$oot_appdir" ]] || die "Application directory '$oot_appdir' as given in specfile $specfile does "\
	"not yet exist. This maybe because you're not in the correct root directory (try cd ..)"\
	" or because it's a nonexisting new application. Don't start new projects off-tree."

export oot_buildname="$2"
export buildprefix="build-"
export oot_outdir="$oot_appdir/$buildprefix$oot_buildname"
mkdir -p $oot_outdir || die "Cannot create $oot_outdir"

export oot_codedir="oot-code" # the subdirectory in $oot_outdir where the files are doubled to
insthash=$(echo $buildscripts | md5sum | fold -w10 | head -n1) # a deterministic hash for this ExaHyPE installation
export qualifiedbuildname="exabuild-${insthash}-$oot_projectname-$oot_buildname"

# Real out of tree support.
if [[ x${OOT_BASE} == "x" ]]; then
	# no base given
	# Ramdisk support: Put all the mirrored files to a fast temporary disk.
	export USE_RAMDISK=${USE_RAMDISK:="Yes"}
	if [[ $USE_RAMDISK == "Yes" ]]; then
		tmproot="/dev/shm"
		export oot_builddir="$tmproot/$qualifiedbuildname"
		mkdir -p $oot_builddir
		[[ -e $oot_outdir/$oot_codedir ]] || ln -s $oot_builddir $oot_outdir/$oot_codedir
		log "Setting up Out-of-Tree copy of all files at ramdisk at $oot_builddir"
	else
		export oot_builddir="$oot_outdir/$oot_codedir"
		mkdir -p $oot_builddir
		log "No OOT_BASE given, setting up the copy of all files at $oot_builddir"
	fi
else
	export oot_builddir="$OOT_BASE/$qualifiedbuildname"
	mkdir -p $oot_builddir
	[[ -e $oot_outdir/$oot_codedir ]] || ln -s $oot_builddir $oot_outdir/$oot_codedir
	log "Using given Out-Of-Tree base $OOT_BASE, so oot_builddir is $oot_builddir"
fi

export logfile="setup.log"

#log "Start logging to $oot_outdir/$logfile. "
#log
#unbuf="stdbuf -i0 -o0 -e0" # turn off buffering in pipe
#exec &> >($unbuf tee "$oot_outdir/$logfile")

# a header for the log
#log "Welcome to the ExaHyPE Out of Tree compiler"
#log

rsync="rsync -r --relative --delete --exclude=.git/ --exclude=.svn/"

# assume that $oot_appdir can be quite populated with a mess
rsync_app_excludelist="--exclude=${buildprefix}*/ --exclude=*.log --exclude=*.vtk --exclude=*.csv --exclude=*.mk --exclude=${makesystem_prefix}*"

log "Copying files..."
verbose $rsync $TOOLKIT $oot_builddir/
verbose $rsync $rsync_app_excludelist $oot_appdir $oot_builddir/
cp $specfile $oot_builddir/$oot_appdir/.. # this is hacky.
verbose $rsync $DEPENDENCIES $oot_builddir/

log "Finished copying to $oot_builddir"
log "Overall oot_builddir size is $(du -hs $oot_builddir | head -n1 | awk '{print $1}') in $(find $oot_builddir | wc -l) files"

log "Exporting configuration:"

env | grep -iE 'oot_' | tee $oot_outdir/oot.env


##### Create scripts in the OutOfTree directory

cat << 'BuildScriptEnd' > $oot_outdir/make.sh
#!/bin/bash
# Mini build script wrapper for ExaHyPE out of tree build system.
# Autogenerated by setup-out-of-tree.sh
# @AUTOGENINFO@

set -e
cd "$(dirname "$0")"
source "oot.env"
cd @OOTCODE@
cd @OOT_APPDIR@

export SKIP_TOOLKIT="${SKIP_TOOLKIT:=No}"
export CLEAN="${CLEAN:=Clean}"

time @COMPILE@

# collect compilation result
verbose cp @OOT_BINARYPATH@ .

BuildScriptEnd

cat << 'DeleteScriptEnd' > $oot_outdir/delete.sh
#!/bin/bash
# Mini script to delete ExaHyPE out of tree build structure.
# Autogenerated by setup-out-of-tree.sh
# @AUTOGENINFO@

set -e
cd "$(dirname "$0")"
codedir=$(readlink -f @OOTCODE@)
rm -r $codedir
rm -rf -- "$(pwd -P)"
DeleteScriptEnd

# replace @variables@ in autogenerated scripts
genscripts="$oot_outdir/make.sh $oot_outdir/delete.sh"
chmod 755 $genscripts
repl() { sed -i "$1" $genscripts; }

repl "s#@AUTOGENINFO@#on $(hostname) at $(date) in $(pwd)#"
repl "s#@OOT_APPDIR@#$oot_appdir#"
repl "s#@COMPILE@#$compile#"
repl "s#@OOT_BINARYPATH@#$oot_codedir/$oot_appdir/$oot_binary#"
repl "s#@OOTCODE@#$oot_codedir#"

log "Created scripts $genscripts in oot builddirectory"

