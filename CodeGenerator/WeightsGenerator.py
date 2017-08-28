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


import Backend
import TemplatingUtils
import Utils #matrix operation and build functions


class WeightsGenerator:
    m_context = {}

    # name of generated output file
    m_filenameRoot = "GaussLegendreQuadrature"
    
    # quadrature nodes and weights mapped onto [0,1]
    m_wGPN       = []
    

    def __init__(self, i_context):
        self.m_context = i_context
        self.m_wGPN, _ = Utils.getGaussLegendre(self.m_context['nDof'])
        

    def generateCode(self):
        l_weightsVector      = Utils.vectorPad(self.m_wGPN, self.m_context['nDofPad'] - self.m_context['nDof'])
        self.m_context['weights1'] = l_weightsVector
        self.m_context['w1Size'] = len(self.m_context['weights1'])
        self.m_context['w1_seq'] = range(self.m_context['w1Size'])
        if(self.m_context['nDim'] == 2):
            # weightsVector is wGPN itself
            l_weightsVector      = Utils.vectorPad(self.m_wGPN, Backend.getPadWidth(len(self.m_wGPN)))
            self.m_context['weights2'] = l_weightsVector
            self.m_context['w2Size'] = len(self.m_context['weights2'])
            self.m_context['w2_seq'] = range(self.m_context['w2Size'])

            # all combinations of two weights, written as an 1D array
            l_weightsVector = [self.m_wGPN[i] * self.m_wGPN[j] for i in range(self.m_context['nDof']) for j in range(self.m_context['nDof'])]
            l_weightsVector = Utils.vectorPad(l_weightsVector, Backend.getPadWidth(len(l_weightsVector)))
            self.m_context['weights3'] = l_weightsVector
            self.m_context['w3Size'] = len(self.m_context['weights3'])
            self.m_context['w3_seq'] = range(self.m_context['w3Size'])

        elif(self.m_context['nDim'] == 3):
            # all combinations of two weights, written as an 1D array
            l_weightsVector = [self.m_wGPN[i] * self.m_wGPN[j] for i in range(self.m_context['nDof']) for j in range(self.m_context['nDof'])]
            l_weightsVector      = Utils.vectorPad(l_weightsVector, Backend.getPadWidth(len(l_weightsVector)))
            self.m_context['weights2'] = l_weightsVector
            self.m_context['w2Size'] = len(self.m_context['weights2'])
            self.m_context['w2_seq'] = range(self.m_context['w2Size'])

            # all combination of three weights, written as an 1D array
            l_weightsVector = [self.m_wGPN[i] * self.m_wGPN[j] * self.m_wGPN[k] for i in range(self.m_context['nDof']) for j in range(self.m_context['nDof']) for k in range(self.m_context['nDof'])]
            l_weightsVector      = Utils.vectorPad(l_weightsVector, Backend.getPadWidth(len(l_weightsVector)))
            self.m_context['weights3'] = l_weightsVector
            self.m_context['w3Size'] = len(self.m_context['weights3'])
            self.m_context['w3_seq'] = range(self.m_context['w3Size'])
            
        else:
            print("WeightsGenerator.__generateWeightsCombinations(): nDim not supported")

        
        self.m_context['maxOrder'] = 9 #inclusive
        self.m_context['GL_seq'] = range(1,self.m_context['maxOrder']+2)
        self.m_context['wGPN'] = [[] for _ in range(self.m_context['maxOrder']+2)]
        self.m_context['xGPN'] = [[] for _ in range(self.m_context['maxOrder']+2)]
        self.m_context['GPN_seq'] = [[] for _ in range(self.m_context['maxOrder']+2)]
        for i in self.m_context['GL_seq']:
            self.m_context['wGPN'][i], self.m_context['xGPN'][i] = Utils.getGaussLegendre(i)
            self.m_context['GPN_seq'][i] = range(i)
        
        #generate files 
        TemplatingUtils.renderAsFile('GaussLegendreQuadrature_h.template',   self.m_filenameRoot+'.h',   self.m_context)
        TemplatingUtils.renderAsFile('GaussLegendreQuadrature_cpp.template', self.m_filenameRoot+'.cpp', self.m_context)

        

        
    
    