import numpy as np
import Backend
import re

class WeightsGenerator:
    # order of the approximation polynomial
    m_order      = -1

    # number of dimensions we simulate
    m_nDim       = -1

    # quadrature nodes mapped onto [0,1]
    m_wGPN       = []

    # SP, DP
    m_precision  = ''

    # name of generated output file
    m_filename = "Weights.h"


    def __init__(self, i_config, i_precision):
        self.m_order     = i_config['nDof']-1
        self.m_nDim      = i_config['nDim']
        self.m_precision = i_precision



        # compute the Gauss-Legendre weights
        _, w = np.polynomial.legendre.leggauss(self.m_order+1)
        # map onto [0,1]
        self.m_wGPN = 0.5*w


    def generateCode(self):
        l_sourceFile = open(self.m_filename, 'a')
        l_includeGuard = '#ifndef _EXAHYPE_KERNELS_ADERDG_OPTIMISED_WEIGHTS_H_\n'   \
                         '#define _EXAHYPE_KERNELS_ADERDG_OPTIMISED_WEIGHTS_H_\n\n'
        l_sourceFile.write(l_includeGuard)
        l_sourceFile.close()

        self.__generateWeightsCombinations()

        # close include guard
        l_sourceFile = open(self.m_filename, 'a')
        l_sourceFile.write('#endif /* _EXAHYPE_KERNELS_ADERDG_OPTIMISED_WEIGHTS_H_ */')
        l_sourceFile.close()


    def __generateWeightsCombinations(self):
        # We need two kinds of combinations
        # (a) aux    = (/ 1.0, wGPN(j), wGPN(k) /)
        #     weight = PRODUCT(aux(1:nDim))
        #     suffix := 2
        # (b) aux    = (/ wGPN(i), wGPN(j), wGPN(k) /)
        #     weight = PRODUCT(aux(1:nDim))
        #     suffix := 3

        if(self.m_nDim == 2):
            # case (a)
            # weightsVector is wGPN itself
            # pad weights vector with zeros
            l_sizeWithoutPadding = np.size(self.m_wGPN) 
            l_padWidth           = Backend.getPadWidth(l_sizeWithoutPadding)
            l_weightsVector      = np.pad(self.m_wGPN, (0.0, l_padWidth), mode='constant')

            self.__writeToFile(l_weightsVector, 2)

            # case (b)
            # all combinations of two weights, written as an 1D array
            l_weightsVector = np.outer(self.m_wGPN, self.m_wGPN).flatten('F')

            # pad this vector with zeros
            l_sizeWithoutPadding = np.size(l_weightsVector)
            l_padWidth           = Backend.getPadWidth(l_sizeWithoutPadding)
            l_weightsVector      = np.pad(l_weightsVector, (0.0, l_padWidth), mode='constant')

            self.__writeToFile(l_weightsVector, 3)

        elif(self.m_nDim == 3):
            # case (a)
            # all combinations of two weights, written as an 1D array
            l_weightsVector = np.outer(self.m_wGPN, self.m_wGPN).flatten('F')

            # pad this vector with zeros
            l_sizeWithoutPadding = np.size(l_weightsVector) 
            l_padWidth           = Backend.getPadWidth(l_sizeWithoutPadding) 
            l_weightsVector      = np.pad(l_weightsVector, (0.0, l_padWidth), mode='constant')

            self.__writeToFile(l_weightsVector, 2)

            # case (b)
            # all combination of three weights, written as an 1D array
            l_weightsVector = np.kron(np.outer(self.m_wGPN, self.m_wGPN), self.m_wGPN).flatten('F')

            # pad this vector with zeros
            l_sizeWithoutPadding = np.size(l_weightsVector)
            l_padWidth           = Backend.getPadWidth(l_sizeWithoutPadding)
            l_weightsVector      = np.pad(l_weightsVector, (0.0, l_padWidth), mode='constant')

            self.__writeToFile(l_weightsVector, 3)

        else:
            print("WeightsGenerator.writeWeightsCombinations(): nDim not supported")


    def __writeToFile(self, i_weightsVector, i_suffix):
        l_vectorSize = np.size(i_weightsVector)

        l_sourceFile = open(self.m_filename, 'a')
 
        l_declaration = "  const DATATYPE kernels::optimised::weights"+str(i_suffix)+"["+str(l_vectorSize)+"] __attribute__((aligned(ALIGNMENT))) = { \n"
        # replace all occurrences of 'DATATYPE' with 'float' and 'double', respectively
        if(self.m_precision=='SP'):
            l_declaration = re.sub(r'\bDATATYPE\b', 'float', l_declaration)
        elif(self.m_precision=='DP'):
            l_declaration = re.sub(r'\bDATATYPE\b', 'double', l_declaration)

        l_content = ""
        for i in range(0, l_vectorSize-1):
            l_content += "    %.16e,\n" % (i_weightsVector[i])
        l_content += "    %.16e\n  }\n\n" % (i_weightsVector[-1])

        l_sourceFile.write(l_declaration)
        l_sourceFile.write(l_content)

        l_sourceFile.close()
