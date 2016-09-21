#!/bin/bash

# This file is part of the ExaHyPE project.
# Copyright (c) 2016  http://exahype.eu
# All rights reserved.
# 
# The project has received funding from the European Union's Horizon 
# 2020 research and innovation programme under grant agreement
# No 671698. For copyrights and licensing, please consult the webpage.
# 
# Released under the BSD 3 Open Source License.
# For the full license text, see LICENSE.txt

MIN_DEGREE=3
MAX_DEGREE=8
CODEGENERATORPATH="../../../../CodeGenerator"

DEGREE=$MIN_DEGREE

echo "Copy kernels and dependencies from Codegenerator"
while [ $DEGREE -le $MAX_DEGREE ]; do
cp $CODEGENERATORPATH/runtime/degree$DEGREE/predictor$DEGREE.cpp ../srcRaw
cp $CODEGENERATORPATH/runtime/degree$DEGREE/GaussLegendreQuadrature$DEGREE.cpp ../srcRaw
cp $CODEGENERATORPATH/runtime/degree$DEGREE/asm_predictor_$DEGREE.c ../srcGen #no need to adapt
let DEGREE=DEGREE+1;
done

echo "Adapt the fetched files using the python script adapt_kernels_for_test.py"
python adapt_kernels_for_test.py