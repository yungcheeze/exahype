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
from jinja2 import Template

class VolumeIntegralGenerator:
    m_config = {}

    # linear/nonlinear
    m_type   = ""

    # name of generated output file
    m_filename = "volumeIntegral.cpp"



    def __init__(self, i_config, i_numerics):
        self.m_config = i_config
        self.m_type   = i_numerics


    def generateCode(self):
        context = {}
        context['nDim'] = self.m_config['nDim']
        context['nVar'] = self.m_config['nVar']
        context['nVarPad'] = Backend.getSizeWithPadding(self.m_config['nVar'])
        context['nDoF'] = self.m_config['nDof']
        context['nDoFPad'] = Backend.getSizeWithPadding(self.m_config['nDof'])
    
        if(self.m_type == 'linear'):
            pass
        else:
            self.__generateNonlinearGemms()
            with open('templates/volumeIntegralNonLinear.template', 'r') as tmp:
                template = Template(tmp.read())
                context['gemm_x'] = 'gemm_'+str(self.m_config['nVar'])+'_'+str(self.m_config['nDof'])+'_'+str(self.m_config['nDof'])+'_lduh_x'
                context['gemm_y'] = 'gemm_'+str(self.m_config['nVar'])+'_'+str(self.m_config['nDof'])+'_'+str(self.m_config['nDof'])+'_lduh_y'
                context['gemm_z'] = 'gemm_'+str(self.m_config['nVar'])+'_'+str(self.m_config['nDof'])+'_'+str(self.m_config['nDof'])+'_lduh_z'
                context['lFhi_padY'] = Backend.getSizeWithPadding(self.m_config['nVar']) \
                                 * Backend.getSizeWithPadding(self.m_config['nDof']**self.m_config['nDim'])
                context['lFhi_padZ'] = 2* Backend.getSizeWithPadding(self.m_config['nVar']) \
                                 * Backend.getSizeWithPadding(self.m_config['nDof']**self.m_config['nDim'])                 
                context['i_seq'] = range(0,self.m_config['nDof'])
                if(self.m_config['nDim'] >= 3):
                    context['j_seq'] = range(0,self.m_config['nDof'])
                else:
                    context['j_seq'] = [0]
                
                with open(self.m_filename, 'w') as out:
                    out.write(template.render(context))


    def __generateNonlinearGemms(self):
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





















