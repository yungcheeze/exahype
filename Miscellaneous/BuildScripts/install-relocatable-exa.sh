#!/bin/bash

# To install the one and only relocatable "exa" tool, use this
# tool. 

# optional command line argument: how to call the binary
exaname="${1:-exa}"

exabin="relocatable-exa.sh"

cd "$(dirname "$0")"

./install-util.sh  "$PWD/$exabin" $exaname || echo "Use '$0 [name of new exa]' to install $exabin with an alias (example: 'exa-disk2')"
