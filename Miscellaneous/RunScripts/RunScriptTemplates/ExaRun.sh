
set -e

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

runenv="run.env"
info "Dumping running environment to $runenv"
env > $runenv

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