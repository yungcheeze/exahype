#!/bin/sh

# to build the optimized kernel for all p orders, you need

export CLEAN="Clean"

# as it distributes code in the exahype/kernels directory.


exa="../../BuildScripts/exa.sh"

$exa polycompile-all OptimisedKernel_Euler

