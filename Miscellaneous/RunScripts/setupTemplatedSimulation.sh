#!/bin/bash
#
# This is a fork of runTemplatedSpecfile.sh which itself was
# a tine wrapper around the exahype binary, handling
# parameters via ENVIRONMENT variables, spec file, simulation
# directory setup and the actual running in one single go,
# moved out and generalized from the ConvergenceStudies.
#
# However it turned out that we more need something like what
# the out-of-tree build system is for the running: A simulation
# setup script which does not do all the heavy lifting on it's
# own but generates a bunch of local scripts instead which can
# be modified or called from different computers, on a queue,
# etc.
#
# The script is, thought, still fully controlled by environmental
# variables and suitable as a replacement for runTemplatedSpecfile.sh.
#
# (c) ExaHyPE, SvenK, 2017

ME="$(basename "$0")"
info () { echo -e $ME: $@; } # print error/info message with name of script
fail () { info $@; exit -1; } # exit after errormessage with name of script
finish () { echo $@; exit 0; } # finish with message happily
silent () { $@ >/dev/null; } # suppress output
verbose() { echo $@; $@; } # show command
newline() { echo; } # semantically improved newline
# put a timestamp before each line, `stdbuf -o0` to avoid 4kB buffering
timestampize() { stdbuf -o0 awk '{ print strftime("%Y-%m-%d %H:%M:%S"), $0; fflush(); }'; }

# directory where simulation are stored in
[[ ! ${SIMDIR} ]] && fail "Need SIMDIR specified to where setup simulation"

# Path to the ExaHyPE executable to use
# possible path is eg. ExaBinary=$(exa root)/$(exa find-binary MHD)
[[ ! ${ExaBinary} ]] && fail "Need ExaBinary as path to the executable"

# path of the spec file to use. This is the template
[[ ! ${ExaSpecfile} ]] && fail "Need ExaSpecfile as path to specification template to use"

# check before touching anything
[[ ! -e "$ExaBinary" ]] && fail "Could not find ExaHyPE binary at '$ExaBinary'."
[[ ! -e "$ExaSpecfile" ]] && fail "Could not find ExaHyPE specificataion file at '$ExaSpecfile'."

# exit script in case of error starting from here
set -e
# enable this is if you run into trouble to see what replacement is done etc.
## set -x

# compose and setup simulation parameter directory
if [ -e "$SIMDIR" ]; then
	info "Wiping existing simulation at '$SIMDIR'"
	rm -r "$SIMDIR";
fi
mkdir -p "$SIMDIR"

# convert possibly relative paths to absolute ones
ExaBinary=$(readlink -f "$ExaBinary")
ExaSpecfile=$(readlink -f "$ExaSpecfile")
SIMDIR=$(readlink -f "$SIMDIR")
LOGFILE=$(readlink -f "$LOGFILE")
# and get the pure filenames
BASE_ExaBinary=$(basename "$ExaBinary")
BASE_ExaSpecfile=$(basename "$ExaSpecfile")

# populate simulation directory with absolute symlinks
ln -s "$ExaBinary" "$SIMDIR/"
cp "$ExaSpecfile" "$SIMDIR/"

cd "$SIMDIR"

if [[ $ExaDetachOutput == "Detach" ]]; then
	info "Calling $BASE_ExaBinary, redirecting output to $LOGFILE"
	exec > >(timestampize > "$LOGFILE") 2>&1
else
	info "Calling $BASE_ExaBinary, logging also to $LOGFILE"
	exec > >(timestampize | tee "$LOGFILE") 2>&1
fi

runenv='ExaRun.env'

# setup a readable, correctly quoted environment file
quoteenv() { perl -e 'foreach $k (sort(keys(%ENV))) { print "$k=\"$ENV{$k}\"\n"; }'; }
quoteenv | grep -iE '^(Exa|Sim)' | tee $runenv


setupenv="setup.env"
info "Dumping setup environment to $setupenv"
env > $setupenv


###### IN THE MIDDLE OF THE WORK MOVING STUFF OVER TO EXARUN.SH


newline
newline

info "Starting ExaHyPE on $(hostname) at $(date):"
verbose $QRUN time ./$BASE_ExaBinary $BASE_ExaSpecfile \
&& info "Finished ExaHyPE successfully" \
|| fail "ExaHyPE binary failed!" 

if silent compgen -G *vtk; then
	info "Packing simulation outcome to results.tar.gzip"
	silent tar cvfz results.tar.gzip *.vtk
else
	info "Not VTK output found, probably a failure in the ExaHyPE run"
fi

finish "Done."

