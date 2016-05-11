import numpy as np
import Backend
import re
import sys
sys.path.insert(0, '../aderdg')
import aderdg

class DGMatrixGenerator:
    m_config = {}

    # order of the approximation polynomial
    m_order      = -1

    # number of degrees of freedom
    m_nDof       = -1

    # number of dimensions we simulate
    m_nDim       = -1

    # quadrature nodes and weights mapped onto [0,1]
    m_xGPN       = []
    m_wGPN       = []

    # SP, DP
    m_precision  = ''

    # linear/nonlinear
    m_numerics  = ''

    # name of generated output file
    m_filename = "DGMatrices.h"


    def __init__(self, i_config, i_precision, i_numerics):
        self.m_order     = i_config['nDof']-1
        self.m_config    = i_config
        self.m_nDim      = i_config['nDim']
        self.m_precision = i_precision
        self.m_numerics    = i_numerics

        # compute the Gauss-Legendre weights
        x, w = np.polynomial.legendre.leggauss(self.m_order+1)
        # map onto [0,1]
        self.m_xGPN = 0.5*(x+1)
        self.m_wGPN = 0.5*w


    def generateCode(self):
        self.__writeIncludeGuard()

        # padding of the system matrices is the same
        l_padHeight = Backend.getPadWidth(self.m_config['nDof'])
        l_padMatrix = np.zeros((l_padHeight, self.m_config['nDof']))

        # Kxi
        Kxi = aderdg.assembleStiffnessMatrix(self.m_xGPN, self.m_wGPN, self.m_order)
        l_linearisedKxi = np.concatenate((Kxi,l_padMatrix),axis=0).flatten('F')
        # Kxi is now
        # [ Kxi ]
        # [0...0]
        # [0...0]
        # stored in column major order
        self.__writeToFile("Kxi", l_linearisedKxi)

        # iK1
        iK1 = np.linalg.inv(aderdg.assembleK1(Kxi, self.m_xGPN, self.m_order))
        l_linearisediK1 = np.concatenate((Kxi,l_padMatrix),axis=0).flatten('F')
        self.__writeToFile("iK1", l_linearisediK1)


        # Cauchy-Kovalevski needs dudx
        if(self.m_numerics == 'linear'):
            MM   = aderdg.assembleMassMatrix(self.m_xGPN, self.m_wGPN, self.m_order)
            dudx = aderdg.assembleDiscreteDerivativeOperator(MM,Kxi)
            # we multiply always from the left -> no transpose
            l_lineariseddudx = np.concatenate((dudx,l_padMatrix),axis=0).flatten('F')
            self.__writeToFile("dudx", l_lineariseddudx)

        # FLCoeff; FRCoeff

        self.__closeIncludeGuard()


    def __writeToFile(self, i_matrixName, i_matrix):
        l_sourceFile = open(self.m_filename, 'a')
        if(self.m_precision=='SP'):
            l_sourceFile.write("const float " + i_matrixName + " __attribute__((aligned(ALIGNMENT))) = {\n")
        else:
            l_sourceFile.write("const double " + i_matrixName + " __attribute__((aligned(ALIGNMENT))) = {\n")
        l_sourceFile.write(',\n'.join(map(str, i_matrix)))
        l_sourceFile.write("};\n\n")
        l_sourceFile.close()


    def __writeIncludeGuard(self):
        l_sourceFile = open(self.m_filename, 'a')
        l_includeGuard = '#ifndef _EXAHYPE_KERNELS_ADERDG_OPTIMISED_DGMATRICES_H_\n'   \
                         '#define _EXAHYPE_KERNELS_ADERDG_OPTIMISED_DGMATRICES_H_\n\n'
        l_sourceFile.write(l_includeGuard)
        l_sourceFile.close()


    def __closeIncludeGuard(self):
        l_sourceFile = open(self.m_filename, 'a')
        l_sourceFile.write('#endif /* _EXAHYPE_KERNELS_ADERDG_OPTIMISED_DGMATRICES_H_ */')
        l_sourceFile.close()