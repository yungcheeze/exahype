#! /usr/bin/env python

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

import re
import string
import os

pathToKernel_raw = '../../../ExaHyPE/kernels'
dir = os.path.dirname(__file__)
pathToKernel = os.path.join(dir,pathToKernel_raw)

print("Fetching and adapting dependencies")

with open(os.path.join(pathToKernel,'GaussLegendreQuadrature.h'), 'r') as inputFile:
    input = inputFile.read(10000000)
    result = input.replace('GAUSSLEGENDRE_H_', 'CODEGENTEST_GEN_GAUSSLEGENDRE_H_')
    
    with open(os.path.join(dir,'Generic/GaussLegendreQuadrature.h'), 'w') as outputFile:
        outputFile.write(result)


with open(os.path.join(pathToKernel,'GaussLegendreQuadrature.cpp'), 'r') as inputFile:
    input = inputFile.read(10000000)
    result = input.replace('#include "kernels/GaussLegendreQuadrature.h"', '#include "GaussLegendreQuadrature.h"')
    
    with open(os.path.join(dir,'Generic/GaussLegendreQuadrature.cpp'), 'w') as outputFile:
        outputFile.write(result)


with open(os.path.join(pathToKernel,'aderdg/optimised/GaussLegendreQuadrature.h'), 'r') as inputFile:
    input = inputFile.read(10000000)
    result = input.replace('_EXAHYPE_KERNELS_ADERDG_OPTIMISED_WEIGHTS_H_', 'CODEGENTEST_GEN_EXAHYPE_KERNELS_ADERDG_OPTIMISED_WEIGHTS_H_')
    
    with open(os.path.join(dir,'Optimised/GaussLegendreQuadrature.h'), 'w') as outputFile:
        outputFile.write(result)

with open(os.path.join(pathToKernel,'aderdg/optimised/GaussLegendreQuadrature.cpp'), 'r') as inputFile:
    input = inputFile.read(10000000)
    result = input.replace('#include "kernels/aderdg/optimised/GaussLegendreQuadrature.h"', '#include "GaussLegendreQuadrature.h"')
    
    with open(os.path.join(dir,'Optimised/GaussLegendreQuadrature.cpp'), 'w') as outputFile:
        outputFile.write(result)


