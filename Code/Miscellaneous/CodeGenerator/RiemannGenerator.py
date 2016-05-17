#!/bin/env python
##
#
#
#--------------------------------------------------------------
# Generates the code for the Riemann solver
# for a specific configuration
#--------------------------------------------------------------
#
#
#
import Backend
from MatmulConfig import MatmulConfig
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

        l_includeStatement = '#include "string.h"\n'                             \
                             '#include "kernels/aderdg/optimised/Kernels.h"\n'   \
                             '#include "kernels/GaussLegendreQuadrature.h"\n\n'

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

        # fix datatype choice
        l_sourceCode=open(self.m_filename).read()
        if('DATATYPE' in l_sourceCode):
            if(self.m_precision=='SP'):
                l_sourceCode=l_sourceCode.replace('DATATYPE', 'float')
            else:
                l_sourceCode=l_sourceCode.replace('DATATYPE', 'double')
            l_file=open(self.m_filename, 'w')
            l_file.write(l_sourceCode)
            l_file.flush()
            l_file.close()


    def __generateAverageStates(self):
        l_file = open(self.m_filename, 'a')

        # uniform length of QavL, QavR
        l_paddedVectorLength = Backend.getSizeWithPadding(self.m_config['nVar'])
        l_file.write('  DATATYPE QavL['+str(l_paddedVectorLength)+'] __attribute__((aligned(ALIGNMENT))) = {0.};\n')
        l_file.write('  DATATYPE QavR['+str(l_paddedVectorLength)+'] __attribute__((aligned(ALIGNMENT))) = {0.};\n\n')

        # we need
        # aux = (/ 1.0, wGPN(j), wGPN(k) /)
        # weight = PRODUCT(aux(1:nDim))
        # correspond to weights2
        #
        # RHS[pressures + padding | momentum_x + padding | momentum_y + padding | momentum_z + padding | energy + padding]
        #                                                                              ^
        #                                                                           sized 0 in 2D
        # Each of the 5 chunks starts at an appropriately aligned address

        l_file.write('#pragma simd\n')
        l_file.write('  for(int i=0;i<'+str(self.m_chunkSize)+';i++) {\n')
        for iVar in range(0, self.m_config['nVar']):
            l_file.write('    QavL['+str(iVar)+'] += kernels::weights2[i] * lQbndL['+str(iVar*self.m_chunkSize)+'+i];\n')
        l_file.write('  }\n\n')

        l_file.write('#pragma simd\n')
        l_file.write('  for(int i=0;i<'+str(self.m_chunkSize)+';i++) {\n')
        for iVar in range(0, self.m_config['nVar']):
            l_file.write('    QavR['+str(iVar)+'] += kernels::weights2[i] * lQbndR['+str(iVar*self.m_chunkSize)+'+i];\n')
        l_file.write('  }\n\n')

        l_file.close()

    def __generateRusanovSolverForNonlinear(self):
        # write #include's and function signature
        self.__writeHeaderForRiemannSolver()

        self.__generateAverageStates()

        l_file = open(self.m_filename, 'a')

        # uniform length of LL, LR
        l_paddedVectorLength = Backend.getSizeWithPadding(self.m_config['nVar'])

        l_file.write('  DATATYPE LL['+str(l_paddedVectorLength)+'] __attribute__((aligned(ALIGNMENT)));\n')
        l_file.write('  PDEEigenvalues(&QavL[0], normalNonZero, &LL[0]);\n')
        l_file.write('  DATATYPE LR['+str(l_paddedVectorLength)+'] __attribute__((aligned(ALIGNMENT)));\n')
        l_file.write('  PDEEigenvalues(&QavR[0], normalNonZero, &LR[0]);\n\n')

        # abs with intrinsics?
        l_file.write('  DATATYPE smax = 0.;\n')
        l_file.write('  for (int ivar = 0; ivar < '+str(self.m_config['nVar'])+'; ivar++) {\n')
        l_file.write('    smax = std::max(smax, std::max(fabs(lambdaL[ivar]), fabs(lambdaR[ivar])));\n')
        l_file.write('  }\n\n')

        # lFbndL(k,j,:) = 0.5 *   ( lFbndR(k,j,:) + lFbndL(k,j,:) ) - 0.5*smax*( lQbndR(k,j,:) - lQbndL(k,j,:) )
        #               = 0.5 *   ( lFbndR(k,j,:) + lFbndL(k,j,:) ) + 0.5*smax*( lQbndL(k,j,:) - lQbndR(k,j,:) )
        #               = 0.5 * [ ( lFbndR(k,j,:) + lFbndL(k,j,:) ) +     smax*( lQbndL(k,j,:) - lQbndR(k,j,:) ) ]
        l_file.write('#pragma simd\n')
        l_file.write('  for(int i=0; i<'+str(self.m_vectorLength)+'; i++) {\n')
        l_file.write('    lQbndL[i] = smax * (lQbndL[i]-lQbndR[i]);\n')
        l_file.write('  }\n')

        l_file.write('#pragma simd\n')
        l_file.write('  for(int i=0; i<'+str(self.m_vectorLength)+'; i++) {\n')
        l_file.write('    lFbndL[i] += lFbndR[i];\n')
        l_file.write('  }\n')

        l_file.write('#pragma simd\n')
        l_file.write('  for(int i=0; i<'+str(self.m_vectorLength)+'; i++) {\n')
        l_file.write('    lFbndL[i] = 0.5 * (lFbndL[i] + lQbndL[]) ;\n')
        l_file.write('  }\n')
        l_file.write('  memcpy(lFbndR,lFbndL,'+str(self.m_vectorLength)+' * sizeof(DATATYPE));\n')


        # write missing closing bracket
        l_file.write('}')
        l_file.close()

    def __generateRusanovSolverForLinear(self):
        # TODO is this still feasible when we have material parameters?
        self.__writeHeaderForRiemannSolver()

        self.__generateAverageStates()

        # TODO
        # average out the material parameters

        l_file = open(self.m_filename, 'a')

        # uniform length of QavL, QavR, LL, LR
        l_paddedVectorLength = Backend.getSizeWithPadding(self.m_config['nVar'])


        # compute averaged state
        l_file.write('  #pragma simd')
        l_file.write('  for(int i=0;i<'+str(l_paddedVectorLength)+';i++) {')
        l_file.write('    QavL[i] += 0.5 * (QavL[i]+QavR[i]);')
        l_file.write('  }')
        # evaluate system matrix
        l_file.write('  // CALL PDEMatrixB(out_Bn, QavL, nv)')

        # write missing closing bracket
        l_file.write('}')
        l_file.close()
