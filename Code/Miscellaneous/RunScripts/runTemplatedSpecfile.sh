#!/bin/bash
#
# This is a tiny wrapper around the exahype binary handling
#  (1) parameters and spec files
#  (2) simulation directory setup
#  (3) actually running
# This was moved and generalized from the ConvergenceStudies.
#
# The script is fully controlled by environmental variables.
#
# -- SK, 2016

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

# How much TBB cores to use for shared memory parallelization?
# we assume a default value here if not given
export ExaTbbCores=${ExaTbbCores:=1}

# Whether to detach from the terminal before executing the binary or not.
# This is useful if you run several instances in parallel.
# Valid values: Detach, Keep
# Default value: Detach
export ExaDetachOutput=${ExaDetachOutput:=Detach}

# Logging all output as soon as simulation was setup:
export LOGFILE=${LOGFILE:="$SIMDIR/run-$(hostname).log"}

# You can use this variable to insert an queuing command like slurms srun
# before launching ExaHyPE.
export QRUN="${QRUN:=}"
# an example for SLURM is:
#   QRUN="srun -n1 --partition=x-men --time=12:00:00"

# removed any bash command line parsing (http://stackoverflow.com/a/14203146)
# in favour of using environment variables only.

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

# setup stuff needed for running the executable
mkdir -p output

# change spec file contents:
function exaspecrepl { sed -i "$1" $BASE_ExaSpecfile; }

# remove any line including @comment@
exaspecrepl '/@comment@/d'
# use this to change other stuff:
#exaspecrepl "s/@FOO@/$BARc/g"
# there are border cases when this does not work, eg. when $ExaFoo
# contains newlines. Then you will get weird errors by sed.

# find and log all "exa" variables, omit typical $PWD and $OLDPWD.
info "Replacing the following parameters in the specfile:"
# loop over all exa variables to replace them in the specfile
exavars=$(compgen -v | grep -i exa)
for var in $exavars; do
	eval value=\$$var;
	replpattern="@${var}@"
	# ignore case sensitivity in @EXABLA@ or @ExaBla@
	exaspecrepl "s|$replpattern|$value|ig" # Delimiter "|" must not be in $value
	# Consistent proof that we have replaced this parameter
	echo "$var=$value"
done | tee parameters.env

# was before a seperate command (error prone):
#env | grep -iE 'exa|sim' | grep -vE "PWD|OLDPWD" > | tee parameters.env

newline
newline

complete="complete.env"
info "Dumping complete environment to $complete"
env > $complete

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

