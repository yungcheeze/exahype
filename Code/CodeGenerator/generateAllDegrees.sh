#!/bin/bash
##
# @file This file is part of the ExaHyPE project.
# @author ExaHyPE Group (exahype@lists.lrz.de)
#
# @section LICENSE
#
# Copyright (c) 2016  http://exahype.eu
# All rights reserved.
#
# The project has received funding from the European Union's Horizon 
# 2020 research and innovation programme under grant agreement
# No 671698. For copyrights and licensing, please consult the webpage.
#
# Released under the BSD 3 Open Source License.
# For the full license text, see LICENSE.txt
#
#
# @section DESCRIPTION
#
# for runtime tests generate the optimised solver kernels
# for multiple degrees of the DG polynomial
#
# @note
# requires to have loaded the following modules
# module load python/3.3_anaconda_nompi
# module load git


# update code generation backend
cd ../../libxsmm
git pull
make generator

# go back code generator's directory
cd -

# delete subdirectories with generated code
rm -r runtime

# For runtime tests I like to change the filename and append the degree to the filename
APPEND_DEGREE=true

MIN_DEGREE=3
MAX_DEGREE=8

DEGREE=$MIN_DEGREE

while [ $DEGREE -le $MAX_DEGREE ]; do
  # execute code generator
  python Driver.py Euler 5 $DEGREE 3 nonlinear hsw ../../libxsmm --precision=DP
  
  # create subdirectory for each degree
  mkdir -p runtime/degree$DEGREE

  # move generated files into subdirectory
  mv ../ExaHyPE/kernels/aderdg/optimised/* runtime/degree$DEGREE

  if [ "$APPEND_DEGREE" = true ]; then
    cd runtime/degree$DEGREE

    # append order to generated files, <filename>.cpp becomes <filename><order>.cpp
    for file in *.{c,cpph,cpp}; do
      if [[ $file == *.c ]]; then
        mv "$file" "${file/.c/_$DEGREE.c}"
      fi
      if [[ $file == *.cpp ]]; then
        mv "$file" "${file/.cpp/$DEGREE.cpp}"
      fi
      if [[ $file == *.cpph ]]; then
        mv "$file" "${file/.cpph/$DEGREE.cpph}"
      fi
    done

    cd ../..
  fi

  let DEGREE=DEGREE+1;
done


