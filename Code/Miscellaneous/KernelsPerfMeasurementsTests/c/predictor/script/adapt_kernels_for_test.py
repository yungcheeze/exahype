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

min_deg = 3
max_deg = 8
for deg in range(min_deg,max_deg+1):
    print('Degree '+str(deg))
    with open('../srcRaw/predictor'+str(deg)+'.cpp', 'r') as inputFile:
        input = inputFile.read(10000000)
        result = input.replace('kernels::aderdg::optimised::predictor', 'predictor'+str(deg)) \
            .replace('#include "kernels/aderdg/optimised/Kernels.h"\n#include "kernels/aderdg/optimised/GaussLegendreQuadrature.h"\n#include "kernels/aderdg/optimised/asm_predictor.c"', '#include "asm_predictor_'+str(deg)+'.c"\n#include "../glue/optimizedKernels.h"')

        with open('../srcGen/predictor'+str(deg)+'.cpp', 'w') as outputFile:
            outputFile.write(result)

    with open('../srcRaw/GaussLegendreQuadrature'+str(deg)+'.cpp', 'r') as inputFile:
        input = inputFile.read(10000000)
        result = input.replace('kernels::aderdg::optimised::', '') \
            .replace('#include "kernels/aderdg/optimised/GaussLegendreQuadrature.h"', '#include "../glue/matrices.h"') \
            .replace('freeGaussLegendreNodesAndWeights(const std::set<int>& orders)', 'freeGaussLegendreNodesAndWeights'+str(deg)+'()') \
            .replace('initGaussLegendreNodesAndWeights(const std::set<int>& orders)', 'initGaussLegendreNodesAndWeights'+str(deg)+'()') \
            .replace('double*', '//double*')

        #remove not needed matrices from file
        lines = result.split('\n')
        output = ''
        for line in lines:
            if 'gaussLegendreNodes' in line:
                continue
            if 'gaussLegendreWeights' in line:
                continue
            if 'weights2' in line:
                continue
            if 'weights3' in line:
                continue
            output += line+'\n'
            
        with open('../srcGen/GaussLegendreQuadrature'+str(deg)+'.cpp', 'w') as outputFile:
            outputFile.write(output)
            
