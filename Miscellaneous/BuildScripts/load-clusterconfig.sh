#!/bin/bash
#
# Just a neat interface to the cluster-configs directory.
# Usage is like without arguments if you have a hostname file.
# In any case, you should source this script. 
#
# The advantage of this script instead of sourcing directly the cluster
# configuration file is that it can be invoked from any directory.
#
# Usage examples:
#
#  source ..../load-clusterconfig.sh  supermuc-mpi
#  source ..../load-clusterconfig.sh  loewe
#  source ..../load-clusterconfig.sh  
#
# Without argument, will look up ~/.hostname
#

if ! [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
	>&2 echo "This script has to be sourced."
fi

has() { type $@ &>/dev/null; } # a way to check if command is available
SCRIPTDIR=$(dirname "${BASH_SOURCE[0]}")

HOST_INFO_FILE="$HOME/.hostname"
if ! [ -z "$1" ]; then
	CLUSTERNAME="$1"
elif [[ -e "$HOST_INFO_FILE" ]]; then
        CLUSTERNAME="$(< $HOST_INFO_FILE)"
else
	>&2 echo "Usage: load-clusterconfig.sh <NameOfCluster>"
	>&2 echo "Or make sure you have a file $HOST_INFO_FILE"
	#exit -1 # this is sourced.
fi

CLUSTERCONFIG_DIR="$SCRIPTDIR/../ClusterConfigs"
CLUSTERCONFIG="${CLUSTERNAME}.cfg"
if [[ -e $CLUSTERCONFIG_DIR/$CLUSTERCONFIG ]]; then
	echo "Loading cluster configuration for $CLUSTERNAME"
	# Files shall be sourced in local directory
	OLDPWD=$PWD
	cd $CLUSTERCONFIG_DIR
	source $CLUSTERCONFIG
	has module && module list
	cd $OLDPWD
else
	echo "Cluster $CLUSTERNAME detected, but no ClusterConfig present."
	# echo "Create one at $CLUSTERCONFIG if you like."
fi


