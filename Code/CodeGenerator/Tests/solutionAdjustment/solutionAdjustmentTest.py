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
from jinja2 import Template


def prepareSolutionAdjustmentTest(dimension, nvar, order, arch, alignment, regenerateOptimisedKernel):

    pathToKernel_raw = '../../../ExaHyPE/kernels'
    pathToCodeGenerator_raw = '../..'

    dir = os.path.dirname(__file__)
    pathToKernel = os.path.join(dir,pathToKernel_raw)
    pathToCodeGenerator = os.path.join(dir,pathToCodeGenerator_raw)

    #generate optimized code
    if(regenerateOptimisedKernel):
        print("Regenerating optimised kernel")

        generate_command = 'python '+os.path.join(pathToCodeGenerator,'Driver.py')+' Euler '+str(nvar)+' '+str(order)+' '+str(dimension)+' nonlinear '+arch+' '+os.path.join(pathToCodeGenerator,'../../libxsmm')+' --precision=DP'
        print("Regenerating optimised kernel with: "+generate_command)
        os.system(generate_command)

        fetch_command = 'python '+os.path.join(dir,'../Dependencies/fetch_dependencies.py')
        print("Fetching dependencies with: "+fetch_command)
        os.system(fetch_command)

    #generate Makefile
    print("Generating Makefile from template")
    with open(os.path.join(dir,'Makefile_template'), 'r') as tmp:
        template = Template(tmp.read())
        context = {}
        context['nVar'] = nvar
        context['order'] = order
        context['alignment'] = alignment
        context['dimension'] = dimension
        if(arch == 'knl'):
            context['arch'] = '-xMIC-AVX512'
        elif(arch == 'hsw'):
            context['arch'] = '-xCORE-AVX2'
        with open(os.path.join(dir,'Makefile'), 'w') as out:
            out.write(template.render(context))


    #fetch and adapt kernels
    print("Fetching and adapting kernels")

    with open(os.path.join(pathToKernel,'aderdg/generic/c/'+str(dimension)+'d/solutionAdjustment.cpph'), 'r') as inputFile:
        input = inputFile.read(10000000)
        result = input.replace('#include "kernels/GaussLegendreQuadrature.h"', '#include "../../Dependencies/Generic/GaussLegendreQuadrature.h"') \
            .replace('#include "tarch/la/Vector.h"', '') \
            .replace('tarch::la::Vector<DIMENSIONS, double>&', 'double*') \
            .replace('kernels::aderdg::generic::c::solutionAdjustment', 'solutionAdjustment_generic') 
        
        with open(os.path.join(dir,'Generic/solutionAdjustment.cpph'), 'w') as outputFile:
            outputFile.write(result)

    with open(os.path.join(pathToKernel,'aderdg/optimised/solutionAdjustment.cpph'), 'r') as inputFile:
        input = inputFile.read(10000000)
        result = input.replace('#include "kernels/aderdg/optimised/GaussLegendreQuadrature.h"', '#include "../../Dependencies/Optimised/GaussLegendreQuadrature.h"') \
            .replace('#include "kernels/aderdg/optimised/Kernels.h"', '') \
            .replace('tarch::la::Vector<DIMENSIONS,double>&', 'double*') \
            .replace('kernels::aderdg::optimised::solutionAdjustment', 'solutionAdjustment_optimised') 
        
        with open(os.path.join(dir,'Optimised/solutionAdjustment.cpph'), 'w') as outputFile:
            outputFile.write(result)
            

    #Done
    print("Done")