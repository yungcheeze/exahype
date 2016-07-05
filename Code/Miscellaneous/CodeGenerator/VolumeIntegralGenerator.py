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
# Generates the code for the volume integral
# for a specific configuration
#

import Backend
from MatmulConfig import MatmulConfig
import FunctionSignatures
import Utils

class VolumeIntegralGenerator:
    m_config = {}

    # linear/nonlinear
    m_type   = ""

    # name of generated output file
    m_filename = "volumeIntegral.cpp"

    # total length
    m_lduhLength = -1


    def __init__(self, i_config, i_numerics):
        self.m_config = i_config
        self.m_type   = i_numerics

        # without padding
        self.m_lduhLength = self.m_config['nVar']*(self.m_config['nDof']**self.m_config['nDim'])


    def generateCode(self):
        self.__writeHeader()

        if(self.m_type == 'linear'):
            self.__generateLinear()
        else:
            self.__generateNonlinear()

        # write missing closing bracket
        l_sourceFile = open(self.m_filename, 'a')
        l_sourceFile.write('}')
        l_sourceFile.close()


    def __writeHeader(self):
        l_description = '// Solve the volume integral \n\n'

        if(self.m_type == 'linear'):
            l_includeStatement = '#include "kernels/aderdg/optimised/Kernels.h"\n\n'
        else:
            l_includeStatement = '#include "kernels/aderdg/optimised/Kernels.h"\n'                    \
                                 '#include "kernels/aderdg/optimised/DGMatrices.h"\n'                 \
                                 '#include "kernels/aderdg/optimised/GaussLegendreQuadrature.h"\n'    \
                                 '#include <cstring>\n'                                               \
                                 '#include "asm_volumeIntegral.c"\n\n'

        l_functionSignature = FunctionSignatures.getVolumeIntegralSignature()+" {\n"

        l_sourceFile = open(self.m_filename, 'a')
        l_sourceFile.write(l_description)
        l_sourceFile.write(l_includeStatement)
        l_sourceFile.write(l_functionSignature)
        l_sourceFile.close()


    def __generateLinear(self):
        pass


    def __generateNonlinear(self):
        l_sourceFile = open(self.m_filename, 'a')

        # Compute
        # (1) lduh(:,:,j,k) = lduh(:,:,j,k) + MATMUL( lFhi_x(:,:,j,k), TRANSPOSE(Kxi) )*PRODUCT(aux(1:nDim))/dx(1)
        # (2) lduh(:,i,:,k) = lduh(:,i,:,k) + MATMUL( lFhi_y(:,:,i,k), TRANSPOSE(Kxi) )*PRODUCT(aux(1:nDim))/dx(2)
        # (3) lduh(:,i,j,:) = lduh(:,i,j,:) + MATMUL( lFhi_z(:,:,i,j), TRANSPOSE(Kxi) )*PRODUCT(aux(1:nDim))/dx(3)

        # define a sequence of matmul configs
        l_matmulList = []


        #-----------------------------
        # implementation file
        #-----------------------------

        # (1) MATMUL( lFhi_x(:,:,j,k), TRANSPOSE(Kxi) )
        l_matmul_x = MatmulConfig(  # M
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
                                    "lduh_x",                                          \
                                    # prefetching
                                    "nopf",                                            \
                                    # type
                                    "gemm")
        l_matmulList.append(l_matmul_x)

        # (2) MATMUL( lFhi_y(:,:,i,k), TRANSPOSE(Kxi) )
        l_matmul_y = MatmulConfig(  # M
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
                                    self.m_config['nVar']*self.m_config['nDof'],       \
                                    # alpha 
                                    1,                                                 \
                                    # beta
                                    1,                                                 \
                                    # alignment A
                                    0,                                                 \
                                    # alignment C
                                    0,                                                 \
                                    # name
                                    "lduh_y",                                          \
                                    # prefetching
                                    "nopf",                                            \
                                    # type
                                    "gemm")
        l_matmulList.append(l_matmul_y)

        # (3) MATMUL( lFhi_z(:,:,i,j), TRANSPOSE(Kxi) )
        l_matmul_z = MatmulConfig(  # M
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
                                    self.m_config['nVar']*(self.m_config['nDof']**2),  \
                                    # alpha 
                                    1,                                                 \
                                    # beta
                                    1,                                                 \
                                    # alignment A
                                    0,                                                 \
                                    # alignment C
                                    0,                                                 \
                                    # name
                                    "lduh_z",                                          \
                                    # prefetching
                                    "nopf",                                            \
                                    # type
                                    "gemm")
        l_matmulList.append(l_matmul_z)

        Backend.generateAssemblerCode("asm_"+self.m_filename, l_matmulList)


        #-----------------------------
        # driver file
        #-----------------------------

        # size of the padded system matrix
        l_nElems = Backend.getSizeWithPadding(self.m_config['nDof']) * self.m_config['nDof']

        # size of tensor slice used in matmul
        l_sliceSize = Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']

        # The function signature specifies lFhi = [lFhi_x | lFhi_y | lFhi_z].
        l_fluxDofSize = Backend.getSizeWithPadding(self.m_config['nVar']) \
                         * Backend.getSizeWithPadding(self.m_config['nDof']**self.m_config['nDim'])
        l_baseAddr_lFhi_x = 0
        l_baseAddr_lFhi_y = 1*l_fluxDofSize
        l_baseAddr_lFhi_z = 2*l_fluxDofSize

        # init lduh
        # Needed since we process the x-, y- and z-direction simulateously and cannot set a
        # beta parameter to zero. Tuning potential for nDim == 2
        l_sourceFile.write("  memset(lduh, 0, sizeof(lduh)*"+str(self.m_lduhLength)+");\n\n")

        l_sourceFile.write("  // Assume equispaced mesh, dx[0] == dx[1] == dx[2]\n")

        if(self.m_config['nDim'] >= 3):
            kmax = self.m_config['nDof']
        else:
            kmax = 1


        # i-th combination of weigths
        i = 0

        # unroll outer loop
        for k in range(0, kmax):
            for j in range(0, self.m_config['nDof']):
                l_dscalCode = Utils.generateDSCAL("kernels::weights2["+str(i)+"]/dx[0]", 
                                                "kernels::Kxi_T", 
                                                "s_m", 
                                                l_nElems)
                l_sourceFile.write(l_dscalCode+"\n")

                # x-direction
                l_startAddr_lduh =   k * self.m_config['nVar']*(self.m_config['nDof']**2) \
                                + j * self.m_config['nVar']*self.m_config['nDof']
                l_startAddr_lFhi_x =   k * Backend.getSizeWithPadding(self.m_config['nVar'])*(self.m_config['nDof']**2) \
                                    + j * Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']
                l_sourceFile.write(
                    "    "+l_matmul_x.baseroutinename
                        +"(&lFhi["+str(l_baseAddr_lFhi_x + l_startAddr_lFhi_x)+"],"\
                        "&kernels::s_m[0],"\
                        "&lduh["+str(l_startAddr_lduh)+"]);\n")

                # y-direction
                l_startAddr_lduh =   k * self.m_config['nVar']*(self.m_config['nDof']**2) \
                                + j * self.m_config['nVar']
                l_startAddr_lFhi_y =   k * Backend.getSizeWithPadding(self.m_config['nVar'])*(self.m_config['nDof']**2) \
                                    + j * Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']
                l_sourceFile.write(
                    "    "+l_matmul_y.baseroutinename
                        +"(&lFhi["+str(l_baseAddr_lFhi_y + l_startAddr_lFhi_y)+"],"\
                        "&kernels::s_m[0],"\
                        "&lduh["+str(l_startAddr_lduh)+"]);\n")

                if(self.m_config['nDim'] >= 3):
                    # z-direction
                    l_startAddr_lduh =   k * self.m_config['nVar']*self.m_config['nDof'] \
                                    + j * self.m_config['nVar']
                    l_startAddr_lFhi_z =   k * Backend.getSizeWithPadding(self.m_config['nVar'])*(self.m_config['nDof']**2) \
                                        + j * Backend.getSizeWithPadding(self.m_config['nVar'])*self.m_config['nDof']
                    l_sourceFile.write("    "+l_matmul_z.baseroutinename
                                            +"(&lFhi["+str(l_baseAddr_lFhi_z + l_startAddr_lFhi_z)+"],"\
                                            "&kernels::s_m[0],"\
                                            "&lduh["+str(l_startAddr_lduh)+"]);\n")

                l_sourceFile.write("\n")

                # next combination of weights
                i = i+1




















