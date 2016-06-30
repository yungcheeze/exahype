#!/bin/bash

# for runtime tests generate the optimised solver kernels
# for multiple degrees of the DG polynomial

# load modules
module load python/3.3_anaconda_nompi
module load git

# update code generation backend
cd ../../../libxsmm
git pull
make generator

# go back code generator's directory
cd -

# delete subdirectories with generated code
rm -r runtime

MIN_DEGREE=3
MAX_DEGREE=7

DEGREE=$MIN_DEGREE

while [ $DEGREE -le $MAX_DEGREE ]; do
  # execute code generator
  python Driver.py Euler 5 $DEGREE 3 nonlinear hsw ../../../libxsmm --precision=DP
  
  # create subdirectory for each degree
  mkdir -p runtime/degree$DEGREE
  
  # move generated files into subdirectory
  mv ../../ExaHyPE/kernels/aderdg/optimised/* runtime/degree$DEGREE
  
  let DEGREE=DEGREE+1;
done


