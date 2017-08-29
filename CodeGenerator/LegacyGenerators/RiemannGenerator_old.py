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
# Generates the code for the Riemann solver
# for a specific configuration
#
# @warning
# * Intel compiler vectorises the nonlinear Riemann solver perfectly. Packed instructions,
#   aligned, wherever possible.
# * Gnu compiler uses suboptimal avx instructions and gives away half the performance. If
#   we want to have a high-performance code with gcc we need an intrinsics generator which
#   prescribes the instructions.
#

import Backend
import FunctionSignatures


class RiemannGenerator:
    m_config = {}

    # linear/nonlinear; later on maybe different types of Riemann solvers
    m_type   = ""

    # name of generated output file
    m_filename = "riemannSolver.cpph"

    # SP, DP
    m_precision = ''

    # total length of work vectors QbndL/QbndR, FbndL/FbndR
    m_vectorLength = -1

    # we have soa-format. Length of one of the arrays.
    m_chunkSize    = -1


    def __init__(self, i_config, i_numerics, i_precision):
        self.m_config    = i_config
        self.m_type      = i_numerics
        self.m_precision = i_precision

        self.m_chunkSize    = Backend.getSizeWithPadding(self.m_config['nDof']**(self.m_config['nDim']-1))
        self.m_vectorLength = self.m_config['nVar'] * self.m_chunkSize


    def __writeHeaderForRiemannSolver(self):
        l_description = '// Solve the Riemann problems \n\n'

        l_includeStatement = '#include <cstring>\n'                             \
                             '#include <cmath>\n'                                \
                             '#include "kernels/aderdg/optimised/Kernels.h"\n'   \
                             '#include "kernels/aderdg/optimised/DGMatrices.h"\n'\
                             '#include "kernels/aderdg/optimised/GaussLegendreQuadrature.h"\n\n'

        l_functionSignature = FunctionSignatures.getRiemannSolverSignature()+" {\n"

        l_sourceFile = open(self.m_filename, 'a')
        l_sourceFile.write(l_description)
        l_sourceFile.write(l_includeStatement)
        l_sourceFile.write(l_functionSignature)
        l_sourceFile.close()


    def generateCode(self):
        if(self.m_type == 'linear'):
            self.__generateRusanovSolverForLinear()
        else:
            self.__generateRusanovSolverForNonlinear()



    def __generateAverageStates(self):
        l_sourceFile = open(self.m_filename, 'a')

        # uniform length of QavL, QavR
        l_paddedVectorLength = Backend.getSizeWithPadding(self.m_config['nVar'])
        l_sourceFile.write('  double QavL['+str(l_paddedVectorLength)+'] __attribute__((aligned(ALIGNMENT))) = {0.};\n')
        l_sourceFile.write('  double QavR['+str(l_paddedVectorLength)+'] __attribute__((aligned(ALIGNMENT))) = {0.};\n\n')

        # we need
        # aux = (/ 1.0, wGPN(j), wGPN(k) /)
        # weight = PRODUCT(aux(1:nDim))
        # correspond to weights2
        #
        # RHS[pressures + padding | momentum_x + padding | momentum_y + padding | momentum_z + padding | energy + padding]
        #                                                                              ^
        #                                                                           sized 0 in 2D
        # Each of the 5 chunks starts at an appropriately aligned address

        for iVar in range(0, self.m_config['nVar']):
            l_sourceFile.write('#pragma simd\n')
            l_sourceFile.write('  for(int i=0;i<'+str(self.m_chunkSize)+';i++) {\n')
            l_sourceFile.write('    QavL['+str(iVar)+'] += weights2[i] * lQbndL['+str(iVar*self.m_chunkSize)+'+i];\n')
            l_sourceFile.write('  }\n\n')

        for iVar in range(0, self.m_config['nVar']):
            l_sourceFile.write('#pragma simd\n')
            l_sourceFile.write('  for(int i=0;i<'+str(self.m_chunkSize)+';i++) {\n')
            l_sourceFile.write('    QavR['+str(iVar)+'] += weights2[i] * lQbndR['+str(iVar*self.m_chunkSize)+'+i];\n')
            l_sourceFile.write('  }\n\n')

        l_sourceFile.close()

    def __generateRusanovSolverForNonlinear(self):
        # write #include's and function signature
        self.__writeHeaderForRiemannSolver()

        l_sourceFile = open(self.m_filename, 'a')

        # gcc and icc specify distinct ways to inform the compiler about guaranteed alignment
        # gcc: double* arr_ = (double*) __builtin_assume_aligned(a, ALIGNMENT);
        # icc: __assume_aligned(a, ALIGNMENT);
        # the default gcc on the cluster exhibits a well-known bug in alignment assumptions
        # => we skip gcc here
        # do not query __GNUC__ - icc also defines this
        l_sourceFile.write('#ifdef __INTEL_COMPILER\n'\
                           '  __assume_aligned(tmp_bnd, ALIGNMENT);\n'\
                           '  __assume_aligned(weights2, ALIGNMENT);\n'
                           '#endif\n')


        self.__generateAverageStates()

        # uniform length of LL, LR
        l_paddedVectorLength = Backend.getSizeWithPadding(self.m_config['nVar'])

        l_sourceFile.write('  double lambdaL['+str(l_paddedVectorLength)+'] __attribute__((aligned(ALIGNMENT)));\n')
        l_sourceFile.write('  PDEEigenvalues(&QavL[0], normalNonZero, &lambdaL[0]);\n')
        l_sourceFile.write('  double lambdaR['+str(l_paddedVectorLength)+'] __attribute__((aligned(ALIGNMENT)));\n')
        l_sourceFile.write('  PDEEigenvalues(&QavR[0], normalNonZero, &lambdaR[0]);\n\n')

        # abs with intrinsics?
        l_sourceFile.write('  double smax = 0.;\n')
        l_sourceFile.write('  for (int ivar = 0; ivar < '+str(self.m_config['nVar'])+'; ivar++) {\n')
        l_sourceFile.write('    smax = std::max(smax, std::max(fabs(lambdaL[ivar]), fabs(lambdaR[ivar])));\n')
        l_sourceFile.write('  }\n\n')

        # lFbndL(k,j,:) = 0.5 *   ( lFbndR(k,j,:) + lFbndL(k,j,:) ) - 0.5*smax*( lQbndR(k,j,:) - lQbndL(k,j,:) )
        #               = 0.5 *   ( lFbndR(k,j,:) + lFbndL(k,j,:) ) + 0.5*smax*( lQbndL(k,j,:) - lQbndR(k,j,:) )
        #               = 0.5 * [ ( lFbndR(k,j,:) + lFbndL(k,j,:) ) +     smax*( lQbndL(k,j,:) - lQbndR(k,j,:) ) ]
        l_sourceFile.write('#pragma simd\n')
        l_sourceFile.write('  for(int i=0; i<'+str(self.m_vectorLength)+'; i++) {\n')
        l_sourceFile.write('    tmp_bnd[i] = smax * (lQbndL[i]-lQbndR[i]);\n')
        l_sourceFile.write('  }\n')

        l_sourceFile.write('#pragma simd\n')
        l_sourceFile.write('  for(int i=0; i<'+str(self.m_vectorLength)+'; i++) {\n')
        l_sourceFile.write('    lFbndL[i] += lFbndR[i];\n')
        l_sourceFile.write('  }\n')

        l_sourceFile.write('#pragma simd\n')
        l_sourceFile.write('  for(int i=0; i<'+str(self.m_vectorLength)+'; i++) {\n')
        l_sourceFile.write('    lFbndL[i] = 0.5 * (lFbndL[i] + tmp_bnd[i]) ;\n')
        l_sourceFile.write('  }\n')
        l_sourceFile.write('  std::memcpy(lFbndR,lFbndL,'+str(self.m_vectorLength)+' * sizeof(double));\n')


        # write missing closing bracket
        l_sourceFile.write('}')
        l_sourceFile.close()

    def __generateRusanovSolverForLinear(self):
        # TODO is this still feasible when we have material parameters?
        self.__writeHeaderForRiemannSolver()

        self.__generateAverageStates()

        # TODO
        # average out the material parameters

        l_sourceFile = open(self.m_filename, 'a')

        # uniform length of QavL, QavR, LL, LR
        l_paddedVectorLength = Backend.getSizeWithPadding(self.m_config['nVar'])


        # compute averaged state
        l_sourceFile.write('  #pragma simd')
        l_sourceFile.write('  for(int i=0;i<'+str(l_paddedVectorLength)+';i++) {')
        l_sourceFile.write('    QavL[i] += 0.5 * (QavL[i]+QavR[i]);')
        l_sourceFile.write('  }')
        # evaluate system matrix
        l_sourceFile.write('  // CALL PDEMatrixB(out_Bn, QavL, nv)')

        # write missing closing bracket
        l_sourceFile.write('}')
        l_sourceFile.close()
