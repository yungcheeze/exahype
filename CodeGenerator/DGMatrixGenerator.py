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
        FLCoeff, _ = np.array(self.BaseFunc1d(0.0, self.m_xGPN, self.m_order)) #is also F0
        FRCoeff, _ = np.array(self.BaseFunc1d(1.0, self.m_xGPN, self.m_order))
        l_paddedFLCoeff = np.pad(FLCoeff, (0, l_padSize), 'constant')
        l_paddedFRCoeff = np.pad(FRCoeff, (0, l_padSize), 'constant')
        self.m_context['FLCoeff'] = l_paddedFLCoeff
        self.m_context['FRCoeff'] = l_paddedFRCoeff
        
        # Matrices are stored in column major order (so the padding should be on the bottom rows)
        # [ Mat ]
        # [0...0]
        # [0...0]
        
        # Kxi
        Kxi = self.assembleStiffnessMatrix(self.m_xGPN, self.m_wGPN, self.m_order)
        self.m_context['Kxi'] = np.pad(Kxi,((0,l_padSize),(0,0)),'constant').flatten('F')

        # Kxi_T
        Kxi_T = np.transpose(Kxi)
        self.m_context['Kxi_T'] = np.pad(Kxi_T,((0,l_padSize),(0,0)),'constant').flatten('F')

        # iK1
        iK1 = np.transpose(np.linalg.inv(self.assembleK1(Kxi, self.m_xGPN, self.m_order)))
        self.m_context['iK1'] = np.pad(iK1,((0,l_padSize),(0,0)),'constant').flatten('F')

        # dudx
        MM   = self.assembleMassMatrix(self.m_xGPN, self.m_wGPN, self.m_order)
        dudx = self.assembleDiscreteDerivativeOperator(MM,Kxi)
        self.m_context['dudx'] = np.pad(dudx,((0,l_padSize),(0,0)),'constant').flatten('F')
        
        # dudx_T
        dudx_T = np.transpose(dudx)
        self.m_context['dudx_T'] = np.pad(dudx_T,((0,l_padSize),(0,0)),'constant').flatten('F')
        
        
        #generate files 
        TemplatingUtils.renderAsFile('DGMatrices_h.template',   self.m_filenameRoot+'.h',   self.m_context)
        TemplatingUtils.renderAsFile('DGMatrices_cpp.template', self.m_filenameRoot+'.cpp', self.m_context)

        
    # Code taken from:    
        # .. module:: aderdg
        # :platform: Unix, Windows, Mac
        # :synopsis: Provides routines to compute ADER-DG basis functions and operators on the unit cube.
        # .. moduleauthor:: Angelika Schwarz <angelika.schwarz@tum.de>
        # :synopsis: Provides routines to compute ADER-DG basis functions and operators on the unit cube.
    
    def BaseFunc1d(self, xi, xin, N):
        """
        Computes the ADER-DG basis functions and their first derivative.
        
        Args:
           xi:
              The reference element point the basis functions are evaluated at.
              Here, xi refers to the greek letter that is often used as a reference element coordinate.
           xin:
              The reference element nodes corresponding to the nodal basis functions.
           N:
              Order of approximation corresponding to N+1 nodal basis functions.
        Returns:
           phi:
              Basis function values.
           phi_xi:
              First derivatives of the basis functions.
        """
        phi    = [1.]*(N+1) 
        phi_xi = [0.]*(N+1)
        for m in range(0,N+1):
            for j in range(0,N+1):
                if j == m:
                    continue 
                phi[m] = phi[m]*(xi-xin[j])/(xin[m]-xin[j])
            for i in range(0,N+1):
                if i == m:
                    continue
                tmp = 1.;
                for j in range(0,N+1):
                    if j == i:
                        continue
                    if j == m:
                        continue
                    tmp = tmp*(xi-xin[j])/(xin[m]-xin[j])
                phi_xi[m] += tmp/(xin[m]-xin[i])
        return phi, phi_xi    

    def assembleStiffnessMatrix(self, xGPN, wGPN, N):
        """
        Computes the (reference) element stiffness matrix for an approximation of
        order N.

        Args:
           xGPN:
              Gauss-Legendre nodes (N nodes).
           wGPN:
              N Gauss-Legendre weights  (N weights).
           N:
              Order of approximation corresponding to N+1 nodal basis functions.
        Returns:
           K_xi:
              The (reference) element stiffness matrix.
        """
        # init matrix with zero
        Kxi = [[0 for _ in range(N+1)] for _ in range(N+1)]
         
        for i in range(0,N+1):
            phi, phi_xi = self.BaseFunc1d(xGPN[i], xGPN, N)
            for k in range(0,N+1):
                for l in range(0,N+1):
                    Kxi[k][l] += wGPN[i]*phi_xi[k]*phi[l] 
            
        return Kxi
    
    def assembleK1(self, Kxi, xGPN, N):
        """
        Computes the difference between the reference element mass operator 
        evaluated at point xi=1.0 and the element stiffness matrix.
        
        Args:
           K_xi:
              The (reference) element stiffness matrix for a approximation of 
              order N.
           xGPN:
              N Gauss-Legendre nodes (N nodes).
           N:
              Order of approximation corresponding to N+1 nodal basis functions.
        Returns:
           K1:
              <unknown>
        """
        phi1, _ = self.BaseFunc1d(1.0, xGPN, N)
        FRm = [[0 for _ in range(N+1)] for _ in range(N+1)]
        
        for k in range(0, N+1):
            for l in range(0, N+1):
                FRm[k][l] = phi1[k]*phi1[l] 
        
        K1 = np.subtract(FRm,Kxi)
        return K1  
        
        
    def assembleMassMatrix(self, xGPN, wGPN, N):
        """
        Computes the (reference) element mass matrix for an approximation of
        order N.

        Args:
           xGPN:
              Gauss-Legendre nodes (N nodes).
           wGPN:
              N Gauss-Legendre weights (N weights).
           N:
              Order of approximation corresponding to N+1 nodal basis functions.
        Returns:
           M_xi:
              The (reference) element mass matrix.
        """
        # init matrix with zeros
        MM = [[0 for _ in range(N+1)] for _ in range(N+1)]
        
        for i in range(0,N+1):
            phi, _ = self.BaseFunc1d(xGPN[i], xGPN, N)
            for k in range(0,N+1):
                for l in range(0,N+1):
                    MM[k][l] += wGPN[i]*phi[k]*phi[l]
          
        return MM
        
        
    def assembleDiscreteDerivativeOperator(self, MM, Kxi):
        """
        Computes some derivative values for debugging purposes.

        Args:
           MM:
              The (reference) element mass matrix for a approximation of 
              order N.
           Kxi:
              The (reference) element stiffness matrix for a approximation of 
              order N.
           
        Returns:
           dudx:
              Derivative values for debugging purposes.
        """
        dudx = np.dot(np.linalg.inv(MM),np.transpose(Kxi))
        return dudx