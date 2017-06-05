#!/bin/bash
#
# This is a lightweight wrapper around exa build-compile suitable for a
# job script. It does not rely on the environment but instead interprets
# its arguments as variables. Usage is like
#
#   batch-compile APP=EulerFlow CLUSTERCONFIG=iboga-gcc-tbb MODE=Release SHAREDMEM=None
#
# and then it just goes and does the job, with heavy outputting, without
# explicitly needing any further environment variables. Of course you can
# use env vars anyway if you like.
#
# SK, 2017-06-03
#

buildscripts="$(dirname "$0")"
set -e # abort on errors
#set -a # export all variables

fail () { >&2 echo $@; exit -1; } # exit after errormessage with name of script

VARS=$($buildscripts/env2source.pl --export $@) && eval $VARS || fail "Failure parsing $@"
[ -z ${APP+x} ] && fail "Need Appname APP"
[ -z ${CLUSTERCONFIG+x} ] && fail "Need Cluster configuration name CLUSTERCONFIG"

# compose name of simulation: ./env2source.pl --keys=,

EXA=${EXA:=$buildscripts/exa.sh}
eval $($EXA clusterconfig $CLUSTERCONFIG)

# compose name of simulation.
# invoke typical stuff. Todo: Use Indirection instead (http://wiki.bash-hackers.org/syntax/pe#indirection)
# so one really can pass like BUILDNAME="%APP-%MODE-%WHATEVER" to get it expanded.
BUILDNAME="$APP-$MODE-$SHAREDMEM-$DISTRIBUTEDMEM"

exec $EXA build-compile $APP $BUILDNAME
