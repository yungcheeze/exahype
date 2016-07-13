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
# Generates the code for the space-time predictor
# for a specific configuration
#

import Backend
from MatmulConfig import MatmulConfig
import FunctionSignatures
import Utils


class SpaceTimePredictorGenerator:
    m_config = {}

    # linear/nonlinear
    m_type   = ""

    # for accessing the quadrature weights and so forth
    m_order  = ""

    # total length
    m_luhLength = -1
    m_rhsLength = -1

    # soa size
    m_luhChunkSize = -1

    # padded aos size
    m_structSize = -1


    def __init__(self, i_config, i_numerics):
        self.m_config     = i_config
        self.m_type       = i_numerics
        self.m_order      = str(self.m_config['nDof']-1)
        self.m_structSize = Backend.getSizeWithPadding(self.m_config['nVar'])

        # without padding of lduh, luh
        self.m_luhChunkSize = (self.m_config['nDof']**self.m_config['nDim'])
        self.m_luhLength    = self.m_config['nVar'] * self.m_luhChunkSize
        # with padding of lduh, luh
        #self.m_luhChunkSize    = Backend.getSizeWithPadding(self.m_config['nDof']**(self.m_config['nDim']-1))
        #self.m_luhLength       = self.m_config['nVar'] * self.m_luhChunkSize
        #
        # alternatively
        #self.m_luhChunkSize    = self.m_config['nDof']**self.m_config['nDim']
        #self.m_luhLength       = Backend.getSizeWithPadding(self.m_luhChunkSize)

        self.m_rhsLength = self.m_config['nVar'] * (self.m_config['nDof']**(self.m_config['nDim']+1))


    def __writeHeaderForPicardLoop(self, i_pathToFile):
        l_description = '// Predictor \n'                                                                 \
                        '// K_{1} \\cdot \\hat{q}^{n+1} + K_{\\xi} \\cdot \\hat{F}^{n} = F_0 \\cdot \\hat{u} \n' \
                        '// computed as \n'                                                               \
                        '// \\hat{q}^{n+1} = K_{1}^{-1}[F_0 \\cdot \\hat{u} - K_{\\xi} \\cdot \\hat{F}^{n}] \n'

        l_includeStatement = '#include "string.h"\n'                                           \
                             '#include "kernels/aderdg/optimised/Kernels.h"\n'                 \
                             '#include "kernels/aderdg/optimised/DGMatrices.h"\n'              \
                             '#include "kernels/aderdg/optimised/GaussLegendreQuadrature.h"\n' \
                             '#include "kernels/aderdg/optimised/asm_picard.c"\n\n'

        l_functionSignature = FunctionSignatures.getPicardLoopSignature(self.m_config['nDim']) + " {\n"

        l_sourceFile = open(i_pathToFile, 'a')
        l_sourceFile.write(l_description)
        l_sourceFile.write(l_includeStatement)
        l_sourceFile.write(l_functionSignature)
        l_sourceFile.close()


    def __writeHeaderForPredictor(self, i_pathToFile):
        l_description = '// Compute the time-averaged space-time polynomials (integration in time) \n\n'

        l_includeStatement = '#include "kernels/aderdg/optimised/Kernels.h"\n'                 \
                             '#include "kernels/aderdg/optimised/GaussLegendreQuadrature.h"\n' \
                             '#include "kernels/aderdg/optimised/asm_predictor.c"\n\n'

        l_functionSignature = FunctionSignatures.getPredictorSignature()+" {\n"

        l_sourceFile = open(i_pathToFile, 'a')
        l_sourceFile.write(l_description)
        l_sourceFile.write(l_includeStatement)
        l_sourceFile.write(l_functionSignature)
        l_sourceFile.close()


    def __writeHeaderForExtrapolator(self, i_pathToFile):
        l_description = '// Compute the boundary-extrapolated values for Q and F*n \n\n'

        l_includeStatement = '#include "kernels/aderdg/optimised/Kernels.h"\n'                 \
                             '#include "kernels/aderdg/optimised/GaussLegendreQuadrature.h"\n' \
                             '#include "kernels/aderdg/optimised/asm_predictor.c"\n'           \
                             '#include "kernels/aderdg/optimised/asm_scatter.c"\n\n'

        l_functionSignature = FunctionSignatures.getExtrapolatorSignature()+" {\n"

        l_sourceFile = open(i_pathToFile, 'a')
        l_sourceFile.write(l_description)
        l_sourceFile.write(l_includeStatement)
        l_sourceFile.write(l_functionSignature)
        l_sourceFile.close()


    def __writeHeaderForCauchyKovalewski(self, i_pathToFile):
        l_includeStatement = '#include "kernels/aderdg/optimised/Kernels.h" \n'   \
                             '#include "kernels/aderdg/optimised/DGMatrices.h"\n' \
                             '#include "string.h"\n'                              \
                             '#include "kernels/aderdg/optimised/asm_cauchyKovalewski.c"\n\n'
        l_functionSignature = FunctionSignatures.getCauchyKovalewskiSignature()+" {\n"

        l_sourceFile = open(i_pathToFile, 'a')
        l_sourceFile.write(l_includeStatement)
        l_sourceFile.write(l_functionSignature)
        l_sourceFile.close()


    def generateCode(self):
        if(self.m_type == 'nonlinear'):
            self.__generatePicardLoop()
            self.__generatePredictor()
            self.__generateExtrapolator()
        else:
            self.__generateCauchyKovalewski()
            self.__generateExtrapolator()


    def __generatePicardLoop(self):
        l_filename = "picard.cpp"

        # write #include's and function signature
        self.__writeHeaderForPicardLoop(l_filename)

        l_sourceFile = open(l_filename, 'a')

        # gcc and icc specify distinct ways to inform the compiler about guaranteed alignment
        # gcc: double* arr_ = (double*) __builtin_assume_aligned(a, ALIGNMENT);
        # icc: __assume_aligned(a, ALIGNMENT);
        # the default gcc on the cluster exhibits a well-known bug in alignment assumptions
        # => we skip gcc here
        # do not query __GNUC__ - icc also defines this
        l_sourceFile.write('#ifdef __INTEL_COMPILER\n'\
                           '  __assume_aligned(kernels::s_m, ALIGNMENT);\n'\
                           '  __assume_aligned(kernels::F0, ALIGNMENT);\n'\
                           '  __assume_aligned(kernels::Kxi, ALIGNMENT);\n'
                           '#endif\n')

        # initialisation of lqh and rhs0:
        # (1) lqh(iVar,l,i,j,k) = luh(iVar,i,j,k)
        # (2) rhs0(iVar,i,j,k,l) = weights3(i,j,k) * F0 * luh(iVar,i,j,k)

        #-----------------------------------------------------------------------------------
        # soa version. later on we probably want to revive this.
        # Note the ordering of luh.
        #-----------------------------------------------------------------------------------
        # (1) lqh(iVar,l,i,j,k) = luh(i,j,k,iVar)
        #nDOF       = str(self.m_config['nDof'])
        #blockWidth = str(self.m_config['nDof'] * self.m_structSize)
        #
        # (1) lqh(iVar,l,i,j,k) = luh(i,j,k,iVar);
        #if(self.m_config['nDim'] == 2):
        #    l_sourceFile.write( "  for(int i=0;i<"+nDOF+";i++) {\n"   \
        #                        "    for(int j=0;j<"+nDOF+";j++) {\n" \
        #                        "       const int lqh_base_addr = ("+str(self.m_config['nDof']**2)+"*j+"
        #                                                            +str(self.m_config['nDof'])   +"*i)*"
        #                                                            +str(self.m_structSize)+";\n" \
        #                        "       const int luh_addr = i + j*"+str(self.m_config['nDof'])+";\n"
        #                       )
        #    for iVar in range(0, self.m_config['nVar']):
        #        l_sourceFile.write("       lqh[lqh_base_addr+"+str(iVar)+"] = luh["+str(iVar*self.m_luhChunkSize)+"+luh_addr];\n")
        #    l_sourceFile.write( "    }\n"\
        #                        "  }\n"
        #                      )
        #
        #if(self.m_config['nDim'] == 3):
        #    l_sourceFile.write( "  for(int i=0;i<"+nDOF+";i++) {\n"     \
        #                        "    for(int j=0;j<"+nDOF+";j++) {\n"   \
        #                        "      for(int k=0;k<"+nDOF+";k++) {\n" \
        #                        "         const int lqh_base_addr = ("+str(self.m_config['nDof']**3)+"*k+"
        #                                                              +str(self.m_config['nDof']**2)+"*j+"
        #                                                              +str(self.m_config['nDof'])   +"*i)*"
        #                                                              +str(self.m_structSize)+";\n" \
        #                        "         const int luh_addr = i + j*"+str(self.m_config['nDof']) +
        #                                                        "+ k*"+str(self.m_config['nDof']**2)+";\n"
        #                       )
        #    for iVar in range(0, self.m_config['nVar']):
        #        l_sourceFile.write("         lqh[lqh_base_addr+"+str(iVar)+"] = luh["+str(iVar*self.m_luhChunkSize)+"+luh_addr];\n")
        #    l_sourceFile.write( "      }\n"
        #                        "    }\n"
        #                        "  }\n"
        #                      )
        # 2D/3D
        #l_sourceFile.write("  //#pragma omp parallel for\n")
        #l_sourceFile.write("  for(int it=0;it<"+str(self.m_config['nDof']**self.m_config['nDim'])+";it++) {\n" \
        #                   "    const int base_addr = it*"+blockWidth+";\n" \
        #                   "    for(int l=1;l<"+nDOF+";l++) {\n" \
        #                   "      memcpy(&lqh[base_addr+l*"+str(self.m_structSize)+"], &lqh[base_addr],"+str(self.m_structSize)+"*sizeof(double));\n"
        #                   "    }\n" \
        #                   "  }\n"
        #                  )
        #
        # (2) rhs0(iVar,i,j,k,l) = weights3(i,j,k) * F0 * luh(i,j,k,iVar)
        # is missing


        #-----------------------------------------------------------------------------------
        # aos version. Probably only a temporary solution.
        # Copied from the generic code base.
        #-----------------------------------------------------------------------------------
        # (1) lqh(iVar,l,i,j,k) = luh(iVar,i,j,k)
        l_nSpaceDof = self.m_config['nDof']**self.m_config['nDim']

        # recall that lqh is padded whereas luh is not
        l_colWidth   = Backend.getSizeWithPadding(self.m_config['nVar']) # size of lqh(:,l,i,j,k)
        l_blockWidth = l_colWidth*self.m_config['nDof']                  # size of lqh(:,:,i,j,k)

        l_sourceFile.write('  for(int ijk=0;ijk<'+str(l_nSpaceDof)+';ijk++) {\n')
        # replicate luh(:,i,j,k) nDOFt times
        for i in range(0, self.m_config['nDof']):
            l_sourceFile.write('    std::memcpy(&lqh[ijk*'+str(l_blockWidth)+'+'+str(i*l_colWidth)+'], '\
                                               '&luh[ijk*'+str(self.m_config['nVar'])+'], '+\
                                                str(self.m_config['nVar'])+'*sizeof(double));\n')
        # close for loop
        l_sourceFile.write('  }\n')


        # (2) rhs0(iVar,i,j,k,l) = weights3(i,j,k) * F0 * luh(iVar,i,j,k)
        rhs_length = self.m_config['nVar']*(self.m_config['nDof']**(self.m_config['nDim']+1))
        l_sourceFile.write('  double* rhs0 = (double*) _mm_malloc('+str(rhs_length)+'*sizeof(double),ALIGNMENT);\n')
        l_sourceFile.write('  double* rhs  = (double*) _mm_malloc('+str(rhs_length)+'*sizeof(double),ALIGNMENT);\n')

        if(self.m_config['nDim']==2):
            l_sourceFile.write('  for(int j=0;j<'+str(self.m_config['nDof'])+';j++) {\n'\
                               '    for(int i=0;i<'+str(self.m_config['nDof'])+';i++) {\n'\
                               '      double w=kernels::weights1[i]*kernels::weights1[j];\n'\
                               '      for(int iVar=0;iVar<'+str(self.m_config['nVar'])+';iVar++) {\n')
            l_sourceFile.write('        int addr = j*'+str(self.m_config['nVar']*self.m_config['nDof'])+\
                                                 '+i*'+str(self.m_config['nVar'])+\
                                                 '+iVar;\n')
            l_dscalCode = Utils.generateDSCAL('(w*luh[addr])',  \
                                                   'kernels::F0',  \
                                                   'kernels::s_m', \
                                                   Backend.getSizeWithPadding(self.m_config['nDof']))
            l_sourceFile.write(Backend.reindentBlock(l_dscalCode,6))
            l_sourceFile.write('\n')
            # write rhs0(iVar,i,j,:) = s_m(:)
            # unroll loop over time dofs; entries to be updated are evenly spread out
            l_baseAddr = self.m_config['nVar']*self.m_config['nDof']**2
            for l in range(0, self.m_config['nDof']):
                l_sourceFile.write('        rhs0['+str(l*l_baseAddr)+'+addr] = s_m['+str(l)+'];\n')

            # close for loops
            l_sourceFile.write('      }\n'\
                               '    }\n'\
                               '  }\n\n')
        if(self.m_config['nDim']==3):
            l_sourceFile.write('  for(int k=0;k<'+str(self.m_config['nDof'])+';k++) {\n'\
                               '    for(int j=0;j<'+str(self.m_config['nDof'])+';j++) {\n'\
                               '      for(int i=0;i<'+str(self.m_config['nDof'])+';i++) {\n'\
                               '        double w=kernels::weights1[i]*kernels::weights1[j]*kernels::weights1[k];\n'\
                               '        for(int iVar=0;iVar<'+str(self.m_config['nVar'])+';iVar++) {\n')
            l_sourceFile.write('          int addr = k*'+str(self.m_config['nVar']*(self.m_config['nDof']**2))+\
                                                   '+j*'+str(self.m_config['nVar']*self.m_config['nDof'])+\
                                                   '+i*'+str(self.m_config['nVar'])+\
                                                   '+iVar;\n')
            l_dscalCode = Utils.generateDSCAL('(w*luh[addr])',  \
                                                   'kernels::F0',  \
                                                   'kernels::s_m', \
                                                   Backend.getSizeWithPadding(self.m_config['nDof']))
            l_sourceFile.write(Backend.reindentBlock(l_dscalCode,8))
            l_sourceFile.write('\n')
            # write rhs0(iVar,i,j,k,:) = s_m(:)
            # unroll loop over time dofs; entries to be updated are evenly spread out
            # use my scatter operator?
            l_baseAddr = self.m_config['nVar']*self.m_config['nDof']**3
            for l in range(0, self.m_config['nDof']):
                l_sourceFile.write('          rhs0['+str(l*l_baseAddr)+'+addr] = s_m['+str(l)+'];\n')

            # close for loops
            l_sourceFile.write('        }\n'\
                               '      }\n'\
                               '    }\n'\
                               '  }\n\n')

        # define a sequence of matmul configs
        l_matmulList = []

        # discrete Picard iterations
        l_sourceFile.write("  for(int iter=0;iter<"+str(self.m_config['nDof'])+";iter++) {\n")

        # distance in terms of memory addresses between flux in x,y,z direction
        fluxOffset = Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']**(self.m_config['nDim']+1)

        # addresses lFh = [lFh_x | lFh_y | lFh_z]
        l_startAddr_lFh_x = 0*fluxOffset
        l_startAddr_lFh_y = 1*fluxOffset
        l_startAddr_lFh_z = 2*fluxOffset

        # in contrast to the Fortran code, the loop over the time dofs has been split up

        # compute the fluxes
        # Loops are being unrolled - temporary solution. Later on the flux function should be vectorised
        if(self.m_config['nDim'] == 3):
            for l in range(0,self.m_config['nDof']):
                for k in range(0,self.m_config['nDof']):
                    for j in range(0,self.m_config['nDof']):
                        for i in range(0,self.m_config['nDof']):
                            lqh_base_addr = \
                                  k*Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']**3\
                                + j*Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']**2\
                                + i*Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']**1\
                                + l*Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']**0
                            lFh_base_addr = \
                                  l*Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']**3\
                                + k*Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']**2\
                                + j*Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']**1\
                                + i*Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']**0
                            # call PDEFlux3d(Q,f,g,h)
                            l_sourceFile.write("  PDEFlux3d(&lqh["+str(lqh_base_addr)+"],"\
                                                           "&lFh["+str(lFh_base_addr)+"],"\
                                                           "&lFh["+str(lFh_base_addr+fluxOffset)+"],"\
                                                           "&lFh["+str(lFh_base_addr+2*fluxOffset)+"]);\n")

        if(self.m_config['nDim'] == 2):
            for l in range(0,self.m_config['nDof']):
                for j in range(0,self.m_config['nDof']):
                    for i in range(0,self.m_config['nDof']):
                        lqh_base_addr = \
                                  j*Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']**2\
                                + i*Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']**1\
                                + l*Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']**0
                        lFh_base_addr = \
                                  l*Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']**2\
                                + j*Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']**1\
                                + i*Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']**0
                        # call PDEFlux2d(Q,f,g)
                        l_sourceFile.write("  PDEFlux2d(&lqh["+str(lqh_base_addr)+"],"\
                                                        "&lFh["+str(lFh_base_addr)+"],"\
                                                        "&lFh["+str(lFh_base_addr+fluxOffset)+"]);\n")


        # Compute the "derivatives" (contributions of the stiffness matrix)
        # (1) rhs(:,:,j,k,l) = rhs0(:,:,j,k,l) - PRODUCT(aux(1:nDim))*dt/dx(1)*MATMUL( lFh(:,:,j,k,l,1), Kxi )
        # (2) rhs(:,i,:,k,l) = rhs(:,i,:,k,l) - PRODUCT(aux(1:nDim))*dt/dx(2)*MATMUL( lFh(:,i,:,k,l,2), Kxi )
        # (3) rhs(:,i,j,:,l) = rhs(:,i,j,:,l) - PRODUCT(aux(1:nDim))*dt/dx(3)*MATMUL( lFh(:,i,j,:,l,3), Kxi )
        # Incorporate minus sign into dt/dx
        l_sourceFile.write("  // Assume equispaced mesh, dx[0] == dx[1] == dx[2]\n")
        l_sourceFile.write("  const double dtdx = -dt/dx[0];\n\n")

        # (1) rhs(:,:,j,k,l) = rhs0(:,:,j,k,l) + PRODUCT(aux(1:nDim))*dt/dx(1)*MATMUL( lFh(:,:,j,k,l,1), Kxi )
        l_sourceFile.write("  memcpy(rhs,rhs0,"+str(self.m_rhsLength)+"*sizeof(double));\n")
        l_matmul = MatmulConfig(    # M
                                    self.m_config['nVar'],                             \
                                    # N
                                    self.m_config['nDof'],                             \
                                    # K
                                    self.m_config['nDof'],                             \
                                    # LDA
                                    Backend.getSizeWithPadding(self.m_config['nVar']), \
                                    # LDB
                                    Backend.getSizeWithPadding(self.m_config['nDof']), \
                                    # LDC
                                    self.m_config['nVar'],                             \
                                    # alpha
                                    1,                                                 \
                                    # beta
                                    1,                                                 \
                                    # alignment A
                                    0,                                                 \
                                    # alignment C
                                    0,                                                 \
                                    # name
                                    "rhs_x",                                           \
                                    # prefetching
                                    "nopf",                                            \
                                    # type
                                    "gemm")
        l_matmulList.append(l_matmul)

        # write the function calls to the driver file
        l_sourceFile.write("  for(int i=0;i<"+str(self.m_config['nDof']**self.m_config['nDim'])+";i++) {\n")
        l_sourceFile.write(Utils.generateDSCAL("dtdx*kernels::weights3[i]",
                                               "kernels::Kxi",
                                               "s_m", self.m_config['nDof']*Backend.getSizeWithPadding(self.m_config['nDof'])))
        l_matrixSize = self.m_config['nVar']*self.m_config['nDof']
        l_paddedMatrixSize = Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']
        l_sourceFile.write("  "+l_matmul.baseroutinename
                               +"(&lFh["+str(l_startAddr_lFh_x)+"+i*"+str(l_paddedMatrixSize)+"]," \
                                " &kernels::s_m[0],"  \
                                " &rhs[i*"+str(l_matrixSize)+"]);\n\n")
        # close for loop
        l_sourceFile.write("  }\n")


        # (2) rhs(:,i,:,k,l) = rhs(:,i,:,k,l) + PRODUCT(aux(1:nDim))*dt/dx(2)*MATMUL( lFh(:,i,:,k,l,2), Kxi )
        l_matmul = MatmulConfig(    # M
                                    self.m_config['nVar'],                             \
                                    # N
                                    self.m_config['nDof'],                             \
                                    # K
                                    self.m_config['nDof'],                             \
                                    # LDA
                                    Backend.getSizeWithPadding(self.m_config['nVar'])
                                    * self.m_config['nDof'],     \
                                    # LDB
                                    Backend.getSizeWithPadding(self.m_config['nDof']), \
                                    # LDC
                                    self.m_config['nVar'] * self.m_config['nDof'],     \
                                    # alpha
                                    1,                                                 \
                                    # beta
                                    1,                                                 \
                                    # alignment A
                                    0,                                                 \
                                    # alignment C
                                    0,                                                 \
                                    # name
                                    "rhs_y",                                           \
                                    # prefetching
                                    "nopf",                                            \
                                    # type
                                    "gemm")
        l_matmulList.append(l_matmul)

        # write the function calls to the driver file
        # merge k and l into one loop
        l_sourceFile.write("  for(int i=0;i<"+str(self.m_config['nDof']**(self.m_config['nDim']-1))+";i++) {\n")
        l_matrixSize = self.m_config['nVar']*self.m_config['nDof']**2
        l_paddedMatrixSize = Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']**2
        # unroll inner loop (i -> nDOFx)
        for i in range(0, self.m_config['nDof']):
            l_sourceFile.write(Utils.generateDSCAL("dtdx*kernels::weights2[i]*kernels::weights1["+str(i)+"]",
                                                   "kernels::Kxi",
                                                   "s_m", self.m_config['nDof']*Backend.getSizeWithPadding(self.m_config['nDof'])))
            l_sourceFile.write("    "+l_matmul.baseroutinename
                                     +"(&lFh["+str(l_startAddr_lFh_y)+"+"+str(i*Backend.getSizeWithPadding(self.m_config['nVar']))+"+i*"+str(l_paddedMatrixSize)+"],"\
                                      " &kernels::s_m[0],"\
                                      " &rhs["+str(i*self.m_config['nVar'])+"+i*"+str(l_matrixSize)+"]);\n")
        # close for loop
        l_sourceFile.write("  }\n")


        # (3) rhs(:,i,j,:,l) = rhs(:,i,j,:,l) + PRODUCT(aux(1:nDim))*dt/dx(3)*MATMUL( lFh(:,i,j,:,l,3), Kxi )
        l_matmul = MatmulConfig(    # M
                                    self.m_config['nVar'],                             \
                                    # N
                                    self.m_config['nDof'],                             \
                                    # K
                                    self.m_config['nDof'],                             \
                                    # LDA
                                    Backend.getSizeWithPadding(self.m_config['nVar'])
                                    * self.m_config['nDof']**2,     \
                                    # LDB
                                    Backend.getSizeWithPadding(self.m_config['nDof']), \
                                    # LDC
                                    self.m_config['nVar'] * self.m_config['nDof']**2,  \
                                    # alpha
                                    1,                                                 \
                                    # beta
                                    1,                                                 \
                                    # alignment A
                                    0,                                                 \
                                    # alignment C
                                    0,                                                 \
                                    # name
                                    "rhs_z",                                           \
                                    # prefetching
                                    "nopf",                                            \
                                    # type
                                    "gemm")

        if(self.m_config['nDim']>=3):
            l_matmulList.append(l_matmul)

            # write the function calls to the driver file
            l_matrixSize = self.m_config['nVar']*self.m_config['nDof']**3
            l_paddedMatrixSize = Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']**3
            # merge i and j into one loop
            l_sourceFile.write("for(int i=0;i<"+str(self.m_config['nDof']**2)+";i++) {\n")
            # unroll outer loop (l -> nDOFt)
            for l in range(0, self.m_config['nDof']):
                l_sourceFile.write(Utils.generateDSCAL("dtdx*kernels::weights2[i]*kernels::weights1["+str(l)+"]",
                                                       "kernels::Kxi",
                                                       "s_m", self.m_config['nDof']*Backend.getSizeWithPadding(self.m_config['nDof'])))
                l_sourceFile.write("    "+l_matmul.baseroutinename
                                        +"(&lFh["+str(l_startAddr_lFh_z)+"+i*"+str(Backend.getSizeWithPadding(self.m_config['nVar']))+"+"+str(l*l_paddedMatrixSize)+"],"\
                                         " &kernels::s_m[0],"\
                                         " &rhs[i*"+str(self.m_config['nVar'])+"+"+str(l*l_matrixSize)+"]);\n")
            # close for loop
            l_sourceFile.write("  }\n")


        # loop over time DOFs done
        # lqh(:,:,i,j,k) = MATMUL( rhs(:,i,j,k,:), TRANSPOSE(iK1) )
        l_matmul = MatmulConfig(    # M
                                    self.m_config['nVar'],                             \
                                    # N
                                    self.m_config['nDof'],                             \
                                    # K
                                    self.m_config['nDof'],                             \
                                    # LDA
                                    self.m_config['nVar']*(self.m_config['nDof']**self.m_config['nDim']), \
                                    # LDB
                                    Backend.getSizeWithPadding(self.m_config['nDof']), \
                                    # LDC
                                    Backend.getSizeWithPadding(self.m_config['nVar']), \
                                    # alpha
                                    1,                                                 \
                                    # beta
                                    0,                                                 \
                                    # alignment A
                                    0,                                                 \
                                    # alignment C
                                    0,                                                 \
                                    # name
                                    "lqh",                                             \
                                    # prefetching
                                    "nopf",                                            \
                                    # type
                                    "gemm")
        l_matmulList.append(l_matmul)

        # write the function call to the driver file
        # note that the DGmatrices.cpp already stores the transpose of iK1
        for i in range(0, self.m_config['nDof']**self.m_config['nDim']):
            l_sourceFile.write(Utils.generateDSCAL("1./kernels::weights3["+str(i)+"]",
                                                   "kernels::iK1",
                                                   "s_m", self.m_config['nDof']*Backend.getSizeWithPadding(self.m_config['nDof'])))
            l_sourceFile.write("  "+l_matmul.baseroutinename
                                   +"(&rhs["+str(i*self.m_config['nVar'])+"]," \
                                    " &kernels::s_m[0],"  \
                                    " &lqh["+str(i*Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof'])+"]);\n\n")

        Backend.generateAssemblerCode("asm_"+l_filename, l_matmulList)


        # close picard iterations loop
        l_sourceFile.write('  } // for iter\n')

        l_sourceFile.write('  _mm_free(rhs0);\n')
        l_sourceFile.write('  _mm_free(rhs);\n')

        # write missing closing bracket
        l_sourceFile.write('}')

        l_sourceFile.close()


    def __generateCauchyKovalewski(self):
        l_filename = "cauchyKovalewski.cpp"
        l_sourceFile = open(l_filename, 'a')
        # TODO
        l_sourceFile.close()


    def __generatePredictor(self):
        l_filename = "predictor.cpp"

        # write #include's and function signature
        self.__writeHeaderForPredictor(l_filename)

        # define a sequence of matmul configs
        l_matmulList = []

        # let's open the file to which we write our function calls to the assembler code
        l_sourceFile = open(l_filename, 'a')


        # in contrast to the Fortran kernels, split the
        # computation of lqhi and lFhi into two separate loops


        #-----------------------------
        # lqhi
        #-----------------------------

        # (1) lqhi(:,i,j,k) = MATMUL(lqh(:,:,i,j,k), wGPN)
        l_matmul = MatmulConfig(    # M
                                    self.m_config['nVar'],                             \
                                    # N
                                    1,                                                 \
                                    # K
                                    self.m_config['nDof'],                             \
                                    # LDA
                                    Backend.getSizeWithPadding(self.m_config['nVar']), \
                                    # LDB
                                    Backend.getSizeWithPadding(self.m_config['nDof']), \
                                    # LDC
                                    Backend.getSizeWithPadding(self.m_config['nVar']), \
                                    # alpha
                                    1,                                                 \
                                    # beta
                                    0,                                                 \
                                    # alignment A
                                    0,                                                 \
                                    # alignment C
                                    0,                                                 \
                                    # name
                                    "lqhi",                                            \
                                    # prefetching
                                    "nopf",                                            \
                                    # type
                                    "gemv")
        l_matmulList.append(l_matmul)


        # write the driver file
        l_iters = self.m_config['nDof'] ** self.m_config['nDim']

        # variant 1: loop not unrolled
        l_matrixSize = Backend.getSizeWithPadding(self.m_config['nVar']) * self.m_config['nDof']
        l_sourceFile.write("  for(int ijk=0;ijk<"+str(l_iters)+";ijk++)\n")
        l_sourceFile.write("    "+l_matmul.baseroutinename\
                          +'(&lqh[ijk*'+str(l_matrixSize)+'],'+  \
                            '&kernels::weights1[0],'+ \
                            '&lqhi[ijk*'+str(Backend.getSizeWithPadding(self.m_config['nVar']))+']);\n')

        # variant 2: loop unrolled
        #l_baseAddrC = 0
        #l_baseAddrA = 0
        #for it in range(0, l_iters):
        #    l_sourceFile.write("  "+l_matmul.baseroutinename+'(&lqh['+str(l_baseAddrA)+'],'+  \
        #                                                            '&kernels::weights1[0],'+ \
        #                                                     '&lqhi['+str(l_baseAddrC)+'])')
        #    l_baseAddrC = l_baseAddrC + Backend.getSizeWithPadding(self.m_config['nVar'])
        #    l_baseAddrA = l_baseAddrA + Backend.getSizeWithPadding(self.m_config['nVar']) * self.m_config['nDof']

        l_sourceFile.write("\n\n")

        #-----------------------------
        # lFhi
        #-----------------------------

        # (2) lFhi_x(:,i,j,k) = MATMUL(lFh(:,i,j,k,:,1), wGPN)
        # (3) lFhi_y(:,j,i,k) = MATMUL(lFh(:,i,j,k,:,2), wGPN)
        # (4) lFhi_z(:,k,i,j) = MATMUL(lFh(:,i,j,k,:,3), wGPN)
        # The function signature specifies lFhi = [lFhi_x | lFhi_y | lFhi_z].
        l_fluxDofSize = Backend.getSizeWithPadding(self.m_config['nVar']) \
                         * Backend.getSizeWithPadding(self.m_config['nDof']**self.m_config['nDim'])
        l_baseAddr_lFhi_x = 0
        l_baseAddr_lFhi_y = 1*l_fluxDofSize
        l_baseAddr_lFhi_z = 2*l_fluxDofSize

        # lFh = [lFh_x | lFh_y | lFh_z]
        lFh_y_startAddr = 1*Backend.getSizeWithPadding(self.m_config['nVar']) \
                           *(self.m_config['nDof']**(self.m_config['nDim']+1))
        lFh_z_startAddr = 2*Backend.getSizeWithPadding(self.m_config['nVar']) \
                           *(self.m_config['nDof']**(self.m_config['nDim']+1))

        # (2)-(4) lFhi_x(:,i,j,k) = MATMUL(lFh(:,i,j,k,:,1), wGPN)
        l_matmul = MatmulConfig(    # M
                                    self.m_config['nVar'],                             \
                                    # N
                                    1,                                                 \
                                    # K
                                    self.m_config['nDof'],                             \
                                    # LDA
                                    Backend.getSizeWithPadding(self.m_config['nVar']) * (self.m_config['nDof']**self.m_config['nDim']), \
                                    # LDB
                                    Backend.getSizeWithPadding(self.m_config['nDof']), \
                                    # LDC
                                    Backend.getSizeWithPadding(self.m_config['nVar']), \
                                    # alpha
                                    1,                                                 \
                                    # beta
                                    0,                                                 \
                                    # alignment A
                                    0,                                                 \
                                    # alignment C
                                    0,                                                 \
                                    # name
                                    "lFhi",                                            \
                                    # prefetching
                                    "nopf",                                            \
                                    # type
                                    "gemv")
        l_matmulList.append(l_matmul)

        # write the driver file
        # (2) lFhi_x(:,i,j,k) = MATMUL(lFh(:,i,j,k,:,1), wGPN)
        l_sourceFile.write("  // lFhi_x(:,i,j,k) = MATMUL(lFh(:,i,j,k,:,1), wGPN)\n")
        l_baseAddrA = 0
        l_baseAddrC = 0
        for it in range(0, l_iters):
            l_sourceFile.write("  "+l_matmul.baseroutinename+'(&lFh['+str(l_baseAddrA)+'],'+  \
                                                             '&kernels::weights1[0],'+\
                                                             '&lFhi['+str(l_baseAddr_lFhi_x+l_baseAddrC)+']);\n')
            l_baseAddrA = l_baseAddrA + Backend.getSizeWithPadding(self.m_config['nVar'])
            l_baseAddrC = l_baseAddrC + Backend.getSizeWithPadding(self.m_config['nVar'])


        # (3) lFhi_y(:,j,i,k) = MATMUL(lFh(:,i,j,k,:,2), wGPN)
        l_sourceFile.write("  // lFhi_y(:,j,i,k) = MATMUL(lFh(:,i,j,k,:,2), wGPN)\n")

        # 2D/3D compatibility
        if(self.m_config['nDim'] == 3):
            kMax = self.m_config['nDof']
        else:
            kMax = 1

        for k in range(0,kMax):
            for j in range(0, self.m_config['nDof']):
                for i in range(0, self.m_config['nDof']):
                    l_baseAddrA =   Backend.getSizeWithPadding(self.m_config['nVar'])*i\
                                  + Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']*j\
                                  + Backend.getSizeWithPadding(self.m_config['nVar'])*(self.m_config['nDof']**2)*k\
                                  + lFh_y_startAddr
                    l_baseAddrC =   Backend.getSizeWithPadding(self.m_config['nVar'])*j\
                                  + Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']*i\
                                  + Backend.getSizeWithPadding(self.m_config['nVar'])*(self.m_config['nDof']**2)*k
                    l_sourceFile.write("  "+l_matmul.baseroutinename\
                                           +'(&lFh['+str(l_baseAddrA)+'],'+  \
                                            '&kernels::weights1[0],'+\
                                            '&lFhi['+str(l_baseAddr_lFhi_y+l_baseAddrC)+']);\n')


        # (4) lFhi_z(:,k,i,j) = MATMUL(lFh(:,i,j,k,:,3), wGPN)
        l_sourceFile.write("  // lFhi_z(:,k,i,j) = MATMUL(lFh(:,i,j,k,:,3), wGPN)\n")
        if(self.m_config['nDim'] >= 3):
            for k in range(0,self.m_config['nDof']):
                for j in range(0, self.m_config['nDof']):
                    for i in range(0, self.m_config['nDof']):
                        l_baseAddrA =   Backend.getSizeWithPadding(self.m_config['nVar'])*i\
                                  + Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']*j\
                                  + Backend.getSizeWithPadding(self.m_config['nVar'])*(self.m_config['nDof']**2)*k\
                                  + lFh_z_startAddr
                        l_baseAddrC =   Backend.getSizeWithPadding(self.m_config['nVar'])*k\
                                  + Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']*j\
                                  + Backend.getSizeWithPadding(self.m_config['nVar'])*(self.m_config['nDof']**2)*i
                        l_sourceFile.write("  "+l_matmul.baseroutinename\
                                           +'(&lFh['+str(l_baseAddrA)+'],'+  \
                                            '&kernels::weights1[0],'+\
                                            '&lFhi['+str(l_baseAddr_lFhi_z+l_baseAddrC)+']);\n')


        # all matmuls have been collected, now launch code generator backend
        Backend.generateAssemblerCode("asm_"+l_filename, l_matmulList)

        # write missing closing bracket
        l_sourceFile.write('}')
        l_sourceFile.close()



    def __generateExtrapolator(self):
        l_filename = "extrapolatedPredictor.cpp"

        # write #include's and function signature
        self.__writeHeaderForExtrapolator(l_filename)

        # define a sequence of matmul configs
        l_matmulList = []

        # The function signature specifies lFhi = [lFhi_x | lFhi_y | lFhi_z].
        l_fluxDofSize = Backend.getSizeWithPadding(self.m_config['nVar']) \
                         * Backend.getSizeWithPadding(self.m_config['nDof']**self.m_config['nDim'])
        l_baseAddr_lFhi_x = 0
        l_baseAddr_lFhi_y = 1*l_fluxDofSize
        l_baseAddr_lFhi_z = 2*l_fluxDofSize

        # structure of lFbnd, lQbnd
        l_chunkSize       = Backend.getSizeWithPadding(self.m_config['nDof']**(self.m_config['nDim']-1))
        l_vectorLength    = self.m_config['nVar'] * l_chunkSize
        l_startAddr_face1 = 0 * l_vectorLength
        l_startAddr_face2 = 1 * l_vectorLength
        l_startAddr_face3 = 2 * l_vectorLength
        l_startAddr_face4 = 3 * l_vectorLength
        l_startAddr_face5 = 4 * l_vectorLength
        l_startAddr_face6 = 5 * l_vectorLength

        # 2D/3D compatibility
        if(self.m_config['nDim'] >= 3):
            kmax = self.m_config['nDof']
        else:
            kmax = 1

        # driver file
        l_sourceFile = open(l_filename, 'a')


        #------------------------------------------
        # lQbnd
        #------------------------------------------

        # ( 1) lQbnd(:,j,k,1) = MATMUL( lqhi(:,:,j,k), FLCoeff )
        # ( 2) lQbnd(:,j,k,2) = MATMUL( lqhi(:,:,j,k), FRCoeff )
        # ( 3) lQbnd(:,i,k,3) = MATMUL( lqhi(:,i,:,k), FLCoeff )
        # ( 4) lQbnd(:,i,k,4) = MATMUL( lqhi(:,i,:,k), FRCoeff )
        # ( 5) lQbnd(:,i,j,5) = MATMUL( lqhi(:,i,j,:), FLCoeff )
        # ( 6) lQbnd(:,i,j,6) = MATMUL( lqhi(:,i,j,:), FRCoeff )


        l_sourceFile.write("  // #pragma this? works on independent data structures!\n")
        # (1),(2) lQbnd(:,j,k,1) = MATMUL( lqhi(:,:,j,k), FLCoeff )
        l_matmul = MatmulConfig(# M
                                self.m_config['nVar'],                              \
                                # N
                                1,                                                  \
                                # K
                                self.m_config['nDof'],                              \
                                # LDA
                                Backend.getSizeWithPadding(self.m_config['nVar']),  \
                                # LDB
                                self.m_config['nDof'],                              \
                                # LDC
                                Backend.getSizeWithPadding(self.m_config['nVar']),  \
                                # alpha
                                1,                                                  \
                                # beta
                                0,                                                  \
                                # alignment A
                                0,                                                  \
                                # alignment C
                                1,                                                  \
                                # name
                                "lQbnd_x",                                          \
                                # prefetching
                                "nopf",                                             \
                                # type
                                "gemv")
        l_matmulList.append(l_matmul)

        l_matrixSize = Backend.getSizeWithPadding(self.m_config['nVar']) * self.m_config['nDof']
        l_sourceFile.write("  for(int jk=0;jk<"+str(self.m_config['nDof']**(self.m_config['nDim']-1))+";jk++)\n")
        l_sourceFile.write("  "+l_matmul.baseroutinename\
                               +'(&lqhi[jk*'+str(l_matrixSize)+'],'+  \
                                ' &kernels::FLCoeff[0],'+\
                                ' &kernels::tmp_bnd[jk*'+str(Backend.getSizeWithPadding(self.m_config['nVar']))+']);\n')
        l_sourceFile.write("  scatter(&kernels::tmp_bnd[0], &lQbnd["+str(l_startAddr_face1)+"]);\n\n")

        l_sourceFile.write("  for(int jk=0;jk<"+str(self.m_config['nDof']**(self.m_config['nDim']-1))+";jk++)\n")
        l_sourceFile.write("  "+l_matmul.baseroutinename\
                               +'(&lqhi[jk*'+str(l_matrixSize)+'],'+  \
                                ' &kernels::FRCoeff[0],'+\
                                ' &kernels::tmp_bnd[jk*'+str(Backend.getSizeWithPadding(self.m_config['nVar']))+']);\n')
        l_sourceFile.write("  scatter(&kernels::tmp_bnd[0], &lQbnd["+str(l_startAddr_face2)+"]);\n\n")

        # (3),(4) lQbnd(:,i,k,3) = MATMUL( lqhi(:,i,:,k), FLCoeff )
        l_matmul = MatmulConfig(# M
                                self.m_config['nVar'],                              \
                                # N
                                1,                                                  \
                                # K
                                self.m_config['nDof'],                              \
                                # LDA
                                Backend.getSizeWithPadding(self.m_config['nVar'])   \
                                    *self.m_config['nDof'],                         \
                                # LDB
                                self.m_config['nDof'],                              \
                                # LDC
                                Backend.getSizeWithPadding(self.m_config['nVar']),  \
                                # alpha
                                1,                                                  \
                                # beta
                                0,                                                  \
                                # alignment A
                                0,                                                  \
                                # alignment C
                                1,                                                  \
                                # name
                                "lQbnd_y",                                          \
                                # prefetching
                                "nopf",                                             \
                                # type
                                "gemv")
        l_matmulList.append(l_matmul)

        for k in range(0, kmax):
            for i in range(0, self.m_config['nDof']):
                l_startAddr_lqhi =   k * Backend.getSizeWithPadding(self.m_config['nVar'])*(self.m_config['nDof']**2)\
                                   + i * Backend.getSizeWithPadding(self.m_config['nVar'])
                l_startAddr_bnd  =   k * Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']\
                                   + i * Backend.getSizeWithPadding(self.m_config['nVar'])
                l_sourceFile.write("  "+l_matmul.baseroutinename\
                     +'(&lqhi['+str(l_startAddr_lqhi)+'],'+  \
                      ' &kernels::FLCoeff[0],'+\
                      ' &kernels::tmp_bnd['+str(l_startAddr_bnd)+']);\n')

        l_sourceFile.write("  scatter(&kernels::tmp_bnd[0], &lQbnd["+str(l_startAddr_face3)+"]);\n\n")

        for k in range(0, kmax):
            for i in range(0, self.m_config['nDof']):
                l_startAddr_lqhi =   k * Backend.getSizeWithPadding(self.m_config['nVar'])*(self.m_config['nDof']**2)\
                                   + i * Backend.getSizeWithPadding(self.m_config['nVar'])
                l_startAddr_bnd  =   k * Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']\
                                   + i * Backend.getSizeWithPadding(self.m_config['nVar'])
                l_sourceFile.write("  "+l_matmul.baseroutinename\
                     +'(&lqhi['+str(l_startAddr_lqhi)+'],'+  \
                      ' &kernels::FRCoeff[0],'+\
                      ' &kernels::tmp_bnd['+str(l_startAddr_bnd)+']);\n')

        l_sourceFile.write("  scatter(&kernels::tmp_bnd[0], &lQbnd["+str(l_startAddr_face4)+"]);\n\n")

        # (5),(6) lQbnd(:,i,j,5) = MATMUL( lqhi(:,i,j,:), FLCoeff )
        l_matmul = MatmulConfig(# M
                                self.m_config['nVar'],                             \
                                # N
                                1,                                                 \
                                # K
                                self.m_config['nDof'],                             \
                                # LDA
                                Backend.getSizeWithPadding(self.m_config['nVar'])  \
                                    *(self.m_config['nDof']**2),                   \
                                # LDB
                                self.m_config['nDof'],                             \
                                # LDC
                                Backend.getSizeWithPadding(self.m_config['nVar']), \
                                # alpha
                                1,                                                 \
                                # beta
                                0,                                                 \
                                # alignment A
                                0,                                                 \
                                # alignment C
                                1,                                                 \
                                # name
                                "lQbnd_z",                                         \
                                # prefetching
                                "nopf",                                            \
                                # type
                                "gemv")
        l_matmulList.append(l_matmul)

        if(self.m_config['nDim'] >= 3):
            l_sourceFile.write("  for(int ij=0;ij<"+str(self.m_config['nDof']**(self.m_config['nDim']-1))+";ij++)\n")
            l_sourceFile.write("  "+l_matmul.baseroutinename\
                                   +'(&lqhi[ij*'+str(Backend.getSizeWithPadding(self.m_config['nVar']))+'],'+  \
                                    ' &kernels::FLCoeff[0],'+\
                                    ' &kernels::tmp_bnd[ij*'+str(Backend.getSizeWithPadding(self.m_config['nVar']))+']);\n')
            l_sourceFile.write("  scatter(&kernels::tmp_bnd[0], &lQbnd["+str(l_startAddr_face5)+"]);\n\n")

            l_sourceFile.write("  for(int ij=0;ij<"+str(self.m_config['nDof']**(self.m_config['nDim']-1))+";ij++)\n")
            l_sourceFile.write("  "+l_matmul.baseroutinename\
                                +'(&lqhi[ij*'+str(Backend.getSizeWithPadding(self.m_config['nVar']))+'],'+  \
                                    ' &kernels::FRCoeff[0],'+\
                                    ' &kernels::tmp_bnd[ij*'+str(Backend.getSizeWithPadding(self.m_config['nVar']))+']);\n')
            l_sourceFile.write("  scatter(&kernels::tmp_bnd[0], &lQbnd["+str(l_startAddr_face6)+"]);\n\n")


        #------------------------------------------
        # lFbnd
        #------------------------------------------

        # (1) lFbnd(:,j,k,1) = MATMUL( lFhi_x(:,:,j,k), FLCoeff )
        # (2) lFbnd(:,j,k,2) = MATMUL( lFhi_x(:,:,j,k), FRCoeff )
        # (3) lFbnd(:,i,k,3) = MATMUL( lFhi_y(:,:,i,k), FLCoeff )
        # (4) lFbnd(:,i,k,4) = MATMUL( lFhi_y(:,:,i,k), FRCoeff )
        # (5) lFbnd(:,i,j,5) = MATMUL( lFhi_z(:,:,i,j), FLCoeff )
        # (6) lFbnd(:,i,j,6) = MATMUL( lFhi_z(:,:,i,j), FRCoeff )


        l_sourceFile.write("\n  // #pragma this? works on independent data structures!\n")
        # (1)-(6) lFbnd(:,j,k,1) = MATMUL( lFhi_x(:,:,j,k), FLCoeff )
        l_matmul = MatmulConfig(# M
                                self.m_config['nVar'],                              \
                                # N
                                1,                                                  \
                                # K
                                self.m_config['nDof'],                              \
                                # LDA
                                Backend.getSizeWithPadding(self.m_config['nVar']),  \
                                # LDB
                                self.m_config['nDof'],                              \
                                # LDC
                                Backend.getSizeWithPadding(self.m_config['nVar']),  \
                                # alpha
                                1,                                                  \
                                # beta
                                0,                                                  \
                                # alignment A
                                0,                                                  \
                                # alignment C
                                1,                                                  \
                                # name
                                "lFbnd",                                            \
                                # prefetching
                                "nopf",                                             \
                                # type
                                "gemv")
        l_matmulList.append(l_matmul)

        l_sourceFile.write("  // x-direction\n")
        l_sourceFile.write("  for(int j=0;j<"+str(self.m_config['nDof']**(self.m_config['nDim']-1))+";j++)\n")
        l_sourceFile.write("    "+l_matmul.baseroutinename\
                                 +"(&lFhi["+str(l_baseAddr_lFhi_x)+"+"+str(l_matrixSize)+"*j],"\
                                  " &kernels::FLCoeff[0],"\
                                  " &kernels::tmp_bnd[j*"+str(Backend.getSizeWithPadding(self.m_config['nVar']))+"]);\n")
        l_sourceFile.write("  scatter(&kernels::tmp_bnd[0], &lFbnd["+str(l_startAddr_face1)+"]);\n")

        l_sourceFile.write("  for(int j=0;j<"+str(self.m_config['nDof']**(self.m_config['nDim']-1))+";j++)\n")
        l_sourceFile.write("    "+l_matmul.baseroutinename\
                                 +"(&lFhi["+str(l_baseAddr_lFhi_x)+"+"+str(l_matrixSize)+"*j],"\
                                  " &kernels::FRCoeff[0],"\
                                  " &kernels::tmp_bnd[j*"+str(Backend.getSizeWithPadding(self.m_config['nVar']))+"]);\n")
        l_sourceFile.write("  scatter(&kernels::tmp_bnd[0], &lFbnd["+str(l_startAddr_face2)+"]);\n\n")

        l_sourceFile.write("  // y-direction\n")
        l_sourceFile.write("  for(int j=0;j<"+str(self.m_config['nDof']**(self.m_config['nDim']-1))+";j++)\n")
        l_sourceFile.write("    "+l_matmul.baseroutinename\
                                 +"(&lFhi["+str(l_baseAddr_lFhi_y)+"+"+str(l_matrixSize)+"*j],"\
                                  " &kernels::FLCoeff[0],"\
                                  " &kernels::tmp_bnd[j*"+str(Backend.getSizeWithPadding(self.m_config['nVar']))+"]);\n")
        l_sourceFile.write("  scatter(&kernels::tmp_bnd[0], &lFbnd["+str(l_startAddr_face3)+"]);\n")

        l_sourceFile.write("  for(int j=0;j<"+str(self.m_config['nDof']**(self.m_config['nDim']-1))+";j++)\n")
        l_sourceFile.write("    "+l_matmul.baseroutinename\
                                 +"(&lFhi["+str(l_baseAddr_lFhi_y)+"+"+str(l_matrixSize)+"*j],"\
                                  " &kernels::FRCoeff[0],"\
                                  " &kernels::tmp_bnd[j*"+str(Backend.getSizeWithPadding(self.m_config['nVar']))+"]);\n")
        l_sourceFile.write("  scatter(&kernels::tmp_bnd[0], &lFbnd["+str(l_startAddr_face4)+"]);\n\n")

        if(self.m_config['nDim'] >= 3):
            l_sourceFile.write("  // z-direction\n")
            l_sourceFile.write("  for(int j=0;j<"+str(self.m_config['nDof']**(self.m_config['nDim']-1))+";j++)\n")
            l_sourceFile.write("    "+l_matmul.baseroutinename\
                                    +"(&lFhi["+str(l_baseAddr_lFhi_z)+"+"+str(l_matrixSize)+"*j],"\
                                    " &kernels::FLCoeff[0],"\
                                    " &kernels::tmp_bnd[j*"+str(Backend.getSizeWithPadding(self.m_config['nVar']))+"]);\n")
            l_sourceFile.write("  scatter(&kernels::tmp_bnd[0], &lFbnd["+str(l_startAddr_face5)+"]);\n")

            l_sourceFile.write("  for(int j=0;j<"+str(self.m_config['nDof']**(self.m_config['nDim']-1))+";j++)\n")
            l_sourceFile.write("    "+l_matmul.baseroutinename\
                                    +"(&lFhi["+str(l_baseAddr_lFhi_z)+"+"+str(l_matrixSize)+"*j],"\
                                    " &kernels::FRCoeff[0],"\
                                    " &kernels::tmp_bnd[j*"+str(Backend.getSizeWithPadding(self.m_config['nVar']))+"]);\n")
            l_sourceFile.write("  scatter(&kernels::tmp_bnd[0], &lFbnd["+str(l_startAddr_face6)+"]);\n\n")


        # all matmuls have been collected, now launch code generator backend
        Backend.generateAssemblerCode("asm_"+l_filename, l_matmulList)

        # write missing closing bracket
        l_sourceFile.write('}\n')
        l_sourceFile.close()





