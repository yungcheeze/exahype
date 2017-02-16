#!/bin/env python
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
# Generates the code for the mapping of the DG polynomial
# onto the [0,1] and calls the user-defined function.
# It's simplistic and just picks the 2D or 3D version
# of the generic kernels.
#

import FunctionSignatures
import Backend

class SolutionAdjustmentGenerator:
    m_config = {}

    # degree of the DG polynomial
    m_order = -1

    # name of generated output file
    m_filename = "solutionAdjustment.cpph"


    def __init__(self, i_config):
        self.m_config = i_config
        self.m_order = i_config['nDof']-1


    def generateCode(self):
        # write includes and function signature
        self.__writeHeader()

        l_sourceFile = open(self.m_filename, 'a')
        l_paddedDim = Backend.getSizeWithPadding(self.m_config['nDim'])
        l_sourceFile.write('  double x['+str(l_paddedDim)+'] __attribute__((aligned(ALIGNMENT)));\n')

        # gcc and icc specify distinct ways to inform the compiler about guaranteed alignment
        # gcc: double* arr_ = (double*) __builtin_assume_aligned(a, ALIGNMENT);
        # icc: __assume_aligned(a, ALIGNMENT);
        # the default gcc on the cluster exhibits a well-known bug in alignment assumptions
        # => we skip gcc here
        # do not query __GNUC__ - icc also defines this
        l_sourceFile.write('#ifdef __INTEL_COMPILER\n'\
                           '  __assume_aligned(kernels::weights1, ALIGNMENT)\n'\
                           '#endif\n')

        # aos format
        if(self.m_config['nDim'] == 2):
            l_sourceFile.write('  for(int ii=0;ii<'+str(self.m_config['nDof'])+';ii++) {\n'\
                               '    const double qr = kernels::aderdg::optimised::gaussLegendreNodes['+str(self.m_order)+'][ii];\n'\
                               '    for(int jj=0;jj<'+str(self.m_config['nDof'])+';jj++) {\n'\
                               '      const double qs = kernels::aderdg::optimised::gaussLegendreNodes['+str(self.m_order)+'][jj];\n'\
                               '      x[0] = center[0] + dx[0] * (qr - 0.5);\n'\
                               '      x[1] = center[1] + dx[1] * (qs - 0.5);\n'\
                               '      const double weight = kernels::aderdg::optimised::weights1[ii] * kernels::aderdg::optimised::weights1[jj];\n'\
                               '      const int startIndex = (ii+'+str(self.m_config['nDof'])+'*jj)*'+str(self.m_config['nVar'])+';\n'\
                               '      PDESolutionAdjustment(x, weight, t, dt, &luh[startIndex]);\n'\
                               '    }\n'\
                               '  }\n')
        elif(self.m_config['nDim'] == 3):
            l_sourceFile.write('  for(int ii=0;ii<'+str(self.m_config['nDof'])+';ii++) {\n'\
                               '    const double qr = kernels::aderdg::optimised::gaussLegendreNodes['+str(self.m_order)+'][ii];\n'\
                               '    for(int jj=0;jj<'+str(self.m_config['nDof'])+';jj++) {\n'\
                               '      const double qs = kernels::aderdg::optimised::gaussLegendreNodes['+str(self.m_order)+'][jj];\n'\
                               '      for(int kk=0;kk<'+str(self.m_config['nDof'])+';kk++) {\n'\
                               '        const double qt = kernels::aderdg::optimised::gaussLegendreNodes['+str(self.m_order)+'][kk];\n'\
                               '        x[0] = center[0] + dx[0] * (qr - 0.5);\n'\
                               '        x[1] = center[1] + dx[1] * (qs - 0.5);\n'\
                               '        x[2] = center[2] + dx[2] * (qt - 0.5);\n'\
                               '        const double weight = kernels::aderdg::optimised::weights1[ii] * kernels::aderdg::optimised::weights1[jj] * kernels::aderdg::optimised::weights1[kk];\n'\
                               '        const int startIndex = (ii+'+str(self.m_config['nDof'])+'*jj+'+str(self.m_config['nDof']**2)+'*kk)*'+str(self.m_config['nVar'])+';\n'\
                               '        PDESolutionAdjustment(x, weight, t, dt, &luh[startIndex]);\n'\
                               '      }\n'\
                               '    }\n'\
                               '  }\n')

        # close function
        l_sourceFile.write('}')
        l_sourceFile.close()


    def __writeHeader(self):
        l_includeStatement = '#include "kernels/aderdg/optimised/Kernels.h"\n' \
                             '#include "kernels/aderdg/optimised/GaussLegendreQuadrature.h"\n\n'

        l_functionSignature = FunctionSignatures.getSolutionAdjustmentSignature()+" {\n"

        l_sourceFile = open(self.m_filename, 'a')
        l_sourceFile.write(l_includeStatement)
        l_sourceFile.write(l_functionSignature)
        l_sourceFile.close()