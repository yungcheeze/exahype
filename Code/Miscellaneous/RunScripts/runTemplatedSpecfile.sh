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

ME="$0"
info () { echo -e $ME: $@; } # print error/info message with name of script
fail () { info $@; exit -1; } # exit after errormessage with name of script
finish () { echo $@; exit 0; } # finish with message happily
silent () { $@ >/dev/null; } # suppress output
verbose() { echo $@; $@; } # show command


# directory where simulation are stored in
[[ ! ${SIMDIR} ]] && fail "Need SIMDIR specified to where setup simulation"

# Path to the ExaHyPE executable to use
# possible path is eg. ExaBinary=$(exa root)/$(exa find-binary MHD)
[[ ! ${ExaBinary} ]] && fail "Need ExaBinary as path to the executable"

# path of the spec file to use. This is the template
[[ ! ${ExaSpecfile} ]] && fail "Need ExaSpecfile as path to specification template to use"

# How much TBB cores to use for shared memory parallelization?
# we assume a default value here if not given
export ExaTbbCores=${ExaTbbCores:=1}

# removed any bash command line parsing (http://stackoverflow.com/a/14203146)
# in favour of using environment variables only.

# exit script in case of error starting from here
set -e

# compose and setup simulation parameter directory
if [ -e "$SIMDIR" ]; then
	info "Wiping existing simulation at '$SIMDIR'"
	rm -r "$SIMDIR";
fi
mkdir -p "$SIMDIR"

# log everything starting from here
# exec > >(awk '{ print strftime("%Y-%m-%d %H:%M:%S"), $0; fflush(); }' |  tee "$SIMDIR/run-$(hostname).log") 2>&1
#
# or instead, ONLY put everything to the log starting from here, surpressing stdout.

# stdbuf -o0 to avoid 4kB buffering
LOGFILE="$SIMDIR/run-$(hostname).log"
echo "Redirecting further output to $LOGFILE and executing ExaHyPE in background"
exec > >(stdbuf -o0 awk '{ print strftime("%Y-%m-%d %H:%M:%S"), $0; fflush(); }' > $LOGFILE) 2>&1

# convert possibly relative paths to absolute ones
ExaBinary=$(readlink -f "$ExaBinary")
ExaSpecfile=$(readlink -f "$ExaSpecfile")
# and get the pure filenames
BASE_ExaBinary=$(basename "$ExaBinary")
BASE_ExaSpecfile=$(basename "$ExaSpecfile")

echo "This is $ME on $(hostname) at $(date)"

# populate simulation directory with absolute symlinks
ln -s "$ExaBinary" "$SIMDIR/"
cp "$ExaSpecfile" "$SIMDIR/"

cd "$SIMDIR"

# setup stuff needed for running the executable
mkdir -p output

# change spec file contents:
function exaspecrepl { sed -i "$1" $BASE_ExaSpecfile; }

# remove any line including @comment@
exaspecrepl '/@comment@/d'
# use this to change other stuff:
#exaspecrepl "s/@WIDTH@/$EXASPEC_WIDTH/g"

# find and log all "exa" variables, omit typical $PWD and $OLDPWD.
echo "Running with the following parameters:"
# loop over all exa variables to replace them in the specfile
exavars=$(compeng -v | grep -i exa)
for var in $exavars; do
	eval value=\$$var;
	replpattern="@${var}@"
	# ignore case sensitivity in @EXABLA@ or @ExaBla@
	exaspecrepl "s/$replpattern/$value/ig"
	# Consistent proof that we have replaced this parameter
	echo "$var=$value"
done | tee parameters.env

# was before a seperate command (error prone):
#env | grep -iE 'exa|sim' | grep -vE "PWD|OLDPWD" > | tee parameters.env

echo

# run it (-> Change echo to verbose)
echo time ./$BASE_ExaBinary $BASE_ExaSpecfile \
   && info "Finished ExaHyPE successfully" \
   || fail "ExaHyPE binary failed!" 

info "Packing simulation outcome to results.tar.gzip"
silent tar cvfz results.tar.gzip *.vtk

finish "Done"




