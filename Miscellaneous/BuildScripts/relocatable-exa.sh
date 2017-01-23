#!/bin/bash

# This is a version of "exa" which deals with multiple installations
# of ExaHyPE, ie. multiple checkouts of the repository on the same
# computer, in a different way.
#
# This version, if put to ~/bin/exa or /usr/local/bin/exa, will invoke
# the exa.sh of the specific installation *dependent on where you
# currently are*, ie. if your working directory is in one installation,
# it will invoke this installation's exa.
#
# This is fundamentally different to the approach of having different
# softlinks, for instance
#
# $ ls -l ~/bin/exa-*
#  exa -> /some/root1/Engine-ExaHyPE/Miscellaneous/BuildScripts/exa.sh
#  exa-iboga -> /some/clusterfs2/ExaHyPE-Engine/Miscellaneous/BuildScripts/exa.sh
#  exa-nonconservative -> /some-branch3/nonconservative/Code/Miscellaneous/BuildScripts/exa.sh
#  exa-workitp -> /some/clusterfs4/Miscellaneous/BuildScripts/exa.sh
#
# Invoking exa that way will always yield in changing this particular
# exa installation while invoking this relocatable-exa.sh will always
# distribute to the particular exa (and report about it). 
#
# This was inspired by simfactory "relocation". 

exa="Miscellaneous/BuildScripts/exa.sh"
curdir="$PWD"

while [[ $PWD != / ]]; do
    cmd="$PWD/$exa"
    if [ -f "$cmd" ]; then
        break
    fi
    unset cmd
    cd ..
done

if [ -z "$cmd" ]; then
    >&2 echo "Could not find $exa, searching up from $curdir"
    exit 1
fi

>&2 echo "Using ExaHyPE at $PWD";
cd "$curdir" # recover curdir

# Forward the call
exec $cmd "$@"

