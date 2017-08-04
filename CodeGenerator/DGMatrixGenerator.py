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

import numpy as np
import re
import sys
import os
dir = os.path.dirname(__file__)
sys.path.insert(1, dir+'/../Miscellaneous/aderdg')
import aderdg

class DGMatrixGenerator:
    m_context = {}

    # name of generated output file
    m_filenameRoot = "DGMatrices"
    
    # quadrature nodes and weights mapped onto [0,1]
    m_xGPN       = []
    m_wGPN       = []
    
    m_order = -1

    def __init__(self, i_context):
        self.m_context = i_context
        
        self.m_order = self.m_context['nDof']-1
        # compute the Gauss-Legendre weights
        x, w = np.polynomial.legendre.leggauss(self.m_context['nDof'])
        # map onto [0,1]
        self.m_xGPN = 0.5*(x+1)
        self.m_wGPN = 0.5*w


    def generateCode(self):
        l_padSize = self.m_context['nDofPad'] - self.m_context['nDof']
        self.m_context['nDofPad_seq'] = range(self.m_context['nDofPad'])
        self.m_context['nDofPadTimesnDof_seq'] = range(self.m_context['nDofPad']*self.m_context['nDof'])

        # [FLCoeff 0...0]; [FRCoeff 0...0];
        # right now FLCoeff, FRCoeff no pad (gives no benefit w.r.t libxsmm)
        FLCoeff, _ = np.array(aderdg.BaseFunc1d(0.0, self.m_xGPN, self.m_order))
        FRCoeff, _ = np.array(aderdg.BaseFunc1d(1.0, self.m_xGPN, self.m_order))
        l_paddedFLCoeff = np.pad(FLCoeff, (0, l_padSize), 'constant')
        l_paddedFRCoeff = np.pad(FRCoeff, (0, l_padSize), 'constant')
        self.m_context['FLCoeff'] = l_paddedFLCoeff
        self.m_context['FRCoeff'] = l_paddedFRCoeff

        # F0 with padding for DSCAL
        F0, _ = np.array(aderdg.BaseFunc1d(0.0, self.m_xGPN, self.m_order))
        l_paddedF0 = np.pad(F0, (0, l_padSize), 'constant')
        self.m_context['F0'] = l_paddedF0
        
        # Matrices are stored in column major order (so the padding should be on the bottom rows)
        # [ Mat ]
        # [0...0]
        # [0...0]
        
        # Kxi
        Kxi = aderdg.assembleStiffnessMatrix(self.m_xGPN, self.m_wGPN, self.m_order)
        self.m_context['Kxi'] = np.pad(Kxi,((0,l_padSize),(0,0)),'constant').flatten('F')

        # Kxi_T
        Kxi_T = np.transpose(Kxi)
        self.m_context['Kxi_T'] = np.pad(Kxi_T,((0,l_padSize),(0,0)),'constant').flatten('F')

        # iK1
        iK1 = np.transpose(np.linalg.inv(aderdg.assembleK1(Kxi, self.m_xGPN, self.m_order)))
        self.m_context['iK1'] = np.pad(iK1,((0,l_padSize),(0,0)),'constant').flatten('F')

        # dudx
        MM   = aderdg.assembleMassMatrix(self.m_xGPN, self.m_wGPN, self.m_order)
        dudx = aderdg.assembleDiscreteDerivativeOperator(MM,Kxi)
        self.m_context['dudx'] = np.pad(dudx,((0,l_padSize),(0,0)),'constant').flatten('F')
        
        # dudx_T
        dudx_T = np.transpose(dudx)
        self.m_context['dudx_T'] = np.pad(dudx_T,((0,l_padSize),(0,0)),'constant').flatten('F')
        
        
        #generate files 
        TemplatingUtils.renderAsFile('DGMatrices_h.template',   self.m_filenameRoot+'.h',   self.m_context)
        TemplatingUtils.renderAsFile('DGMatrices_cpp.template', self.m_filenameRoot+'.cpp', self.m_context)
