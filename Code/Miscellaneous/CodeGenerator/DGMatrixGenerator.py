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

    # linear/nonlinear
    m_numerics  = ''

    # name of generated output file
    m_sourceName = "DGMatrices.cpp"
    m_headerName = "DGMatrices.h"


    def __init__(self, i_config, i_numerics):
        self.m_order     = i_config['nDof']-1
        self.m_config    = i_config
        self.m_nDim      = i_config['nDim']
        self.m_numerics  = i_numerics

        # compute the Gauss-Legendre weights
        x, w = np.polynomial.legendre.leggauss(self.m_order+1)
        # map onto [0,1]
        self.m_xGPN = 0.5*(x+1)
        self.m_wGPN = 0.5*w


    def generateCode(self):
        # padding of the system matrices is the same
        l_padHeight = Backend.getPadWidth(self.m_config['nDof'])
        l_padMatrix = np.zeros((l_padHeight, self.m_config['nDof']))

        l_matrices = {}

        # Kxi
        Kxi = aderdg.assembleStiffnessMatrix(self.m_xGPN, self.m_wGPN, self.m_order)
        l_linearisedKxi = np.concatenate((Kxi,l_padMatrix),axis=0).flatten('F')
        # Kxi is now
        # [ Kxi ]
        # [0...0]
        # [0...0]
        # stored in column major order
        l_matrices['Kxi'] = l_linearisedKxi

        # iK1
        iK1 = np.linalg.inv(aderdg.assembleK1(Kxi, self.m_xGPN, self.m_order))
        l_linearisediK1 = np.concatenate((Kxi,l_padMatrix),axis=0).flatten('F')
        l_matrices['iK1'] = l_linearisediK1

        # dudx (only needed for Cauchy-Kovalevski)
        if(self.m_numerics == 'linear'):
            MM   = aderdg.assembleMassMatrix(self.m_xGPN, self.m_wGPN, self.m_order)
            dudx = aderdg.assembleDiscreteDerivativeOperator(MM,Kxi)
            # we multiply always from the left -> no transpose
            l_lineariseddudx  = np.concatenate((dudx,l_padMatrix),axis=0).flatten('F')
            l_matrices['iK1'] = l_lineariseddudx

        # [FLCoeff 0...0]; [FRCoeff 0...0];
        FLCoeff, _ = np.array(aderdg.BaseFunc1D(0.0, self.m_xGPN, self.m_order))
        FRCoeff, _ = np.array(aderdg.BaseFunc1D(1.0, self.m_xGPN, self.m_order))
        l_paddedFLCoeff = np.pad(FLCoeff, [0, self.getPadWidth(self.m_config['nDof'])], mode='constant')
        l_paddedFRCoeff = np.pad(FRCoeff, [0, self.getPadWidth(self.m_config['nDof'])], mode='constant')
        l_matrices['FLCoeff'] = l_paddedFLCoeff
        l_matrices['FRCoeff'] = l_paddedFRCoeff

        # F0 no padding
        F0, _ = np.array(aderdg.BaseFunc1D(0.0, self.m_xGPN, self.m_order))
        l_matrices['F0'] = F0

        self.__generateHeaderFile()
        self.__writeToFile(l_matrices)

    def __generateHeaderFile(self):
        l_sourceFile = open(self.m_headerName, 'a')
        # copied from Dominic's DGMatrices.h; matrices are linearised
        l_sourceFile.write('#ifndef _EXAHYPE_KERNELS_ADERDG_OPTIMISED_DGMATRICES_H_\n'   \
                           '#define _EXAHYPE_KERNELS_ADERDG_OPTIMISED_DGMATRICES_H_\n\n')
        l_sourceFile.write('namespace kernels { \n'     \
                           'void initDGMatrices();\n\n' \
                           'extern double *Kxi;\n'      \
                           'extern double *iK1;\n'      \
                           'extern double *F0; \n'      \
                           'extern double *FLCoeff;\n'  \
                           'extern double *FRCoeff;\n'  \
                           'extern double ***subOutputMatrix;\n\n')
        l_sourceFile.write('#endif /* _EXAHYPE_KERNELS_ADERDG_OPTIMISED_DGMATRICES_H_ */')
        l_sourceFile.close()


    def __writeToFile(self, i_matrices):
        l_sourceFile = open(self.m_sourceName, 'a')
        l_sourceFile.write('#include "'+self.m_headerName+'"\n' \
                           '#include <mm_malloc.h> //g++\n\n')
        l_sourceFile.write('double* kernels::Kxi;\n' \
                           'double* kernels::iK1;\n' \
                           'double* kernels::F0;\n' \
                           'double* kernels::FLCoeff;\n' \
                           'double* kernels::FRCoeff;\n\n' )
                           #'double*** kernels::subOutputMatrix;\n\n')
        l_sourceFile.write('void kernels::initDGMatrices() {\n')

        for matrix in i_matrices:
            l_elemCount = i_matrices[matrix].size
            l_sourceFile.write(str(matrix)+' = (double *) _mm_malloc(sizeof(double)*'+str(l_elemCount)+', ALIGNMENT);\n')

        for matrix in i_matrices:
            for idx in range(0, i_matrices[matrix].size):
                l_sourceFile.write(str(matrix) + '['+str(idx)+'] = %.12e' % i_matrices[matrix].item(idx)+';\n')

        # apparently subOutputMatrix is nowhere in use
        #subOutputMatrix  = aderdg.assembleSubOutputMatrix(self.m_xGPN, self.m_order, self.m_dim)
        #writeMatrixLookupTableInitToFile(out, "", subOutputMatrix, (order+1)**dim, (order+1)**dim, "subOutputMatrix[%d]" % order)

        l_sourceFile.write("};\n\n")
        l_sourceFile.close()
