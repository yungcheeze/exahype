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
    
    #volumeIntegral optimized kernel
    with open('../srcRaw/volumeIntegral'+str(deg)+'.cpp', 'r') as inputFile:
        input = inputFile.read(10000000)
        result = input.replace('kernels::aderdg::optimised::volumeIntegral', 'volumeIntegral'+str(deg)) \
            .replace('#include "kernels/aderdg/optimised/Kernels.h"\n#include "kernels/aderdg/optimised/DGMatrices.h"\n#include "kernels/aderdg/optimised/GaussLegendreQuadrature.h"\n#include <cstring>\n#include "asm_volumeIntegral.c"\n', '#include "asm_volumeIntegral_'+str(deg)+'.c"\n#include <cstring>\n#include "../glue/optimizedKernels.h"') \
            .replace('const tarch::la::Vector<DIMENSIONS,double> &dx', 'const double* dx')

        with open('../srcGen/volumeIntegral'+str(deg)+'.cpp', 'w') as outputFile:
            outputFile.write(result)

    #DGMatrices, remove useless matrices definitions
    with open('../srcRaw/DGMatrices'+str(deg)+'.cpp', 'r') as inputFile:
        input = inputFile.read(10000000)
        result = input.replace('kernels::aderdg::optimised::', '') \
            .replace('#include "kernels/aderdg/optimised/DGMatrices.h"', '#include "../glue/matrices.h"') \
            .replace('freeDGMatrices(const std::set<int>& orders)', 'freeDGMatrices'+str(deg)+'()') \
            .replace('initDGMatrices(const std::set<int>& orders)', 'initDGMatrices'+str(deg)+'()') \
            .replace('double*', '//double*')

        #remove not needed matrices from file
        lines = result.split('\n')
        output = ''
        for line in lines:
            if 'Kxi' in line and not 'Kxi_T' in line:
                continue
            if 'iK1' in line:
                continue
            if 'dudx' in line:
                continue
            if 's_v' in line:
                continue
            if 'tmp_bnd' in line:
                continue
            if 'F0' in line:
                continue
            if 'FLCoeff' in line:
                continue
            if 'FRCoeff' in line:
                continue
            if 'equidistantGridProjector1d' in line:
                continue
            if 'fineGridProjector1d' in line:
                continue
            output += line+'\n'

        with open('../srcGen/DGMatrices'+str(deg)+'.cpp', 'w') as outputFile:
            outputFile.write(output)
    
    #GaussLegendreQuadrature, remove useless matrices definitions
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
            if 'weights1' in line:
                continue
            if 'weights3' in line:
                continue
            output += line+'\n'
            
        with open('../srcGen/GaussLegendreQuadrature'+str(deg)+'.cpp', 'w') as outputFile:
            outputFile.write(output)
            
