import numpy as np
import Backend
import re

class WeightsGenerator:
    m_order      = -1
    m_nDim       = -1
    m_wGPN       = []
    m_precision  = ''
    
    def __init__(self, i_config, i_precision):
        self.m_order     = i_config['nDof']-1
        self.m_nDim      = i_config['nDim']
        self.m_precision = i_precision
        
        # compute the Gauss-Legendre weights
        _, w = np.polynomial.legendre.leggauss(self.m_order+1)
        # map onto [0,1]
        self.m_wGPN = 0.5*w        


    def writeWeightsCombinations(self, i_pathToOutputFile):
        if(self.m_nDim == 2):
            pass
            
        elif(self.m_nDim == 3):
            # all combinations of two weights, written as as 1D array
            l_weightsVector = np.outer(self.m_wGPN, self.m_wGPN).flatten('F')

            # pad this vector with zeros
            l_sizeWithoutPadding = np.size(l_weightsVector) 
            l_padWidth           = Backend.getPadWidth(l_sizeWithoutPadding) 
            l_weightsVector      = np.pad(l_weightsVector, (0,l_padWidth), mode='constant')

            self.__writeToFile(i_pathToOutputFile, l_weightsVector, 2)
            

            print(l_weightsVector)
        else:
            print("WeightsGenerator.generatePredictor(): nDim not supported")
            

    def __writeToFile(self, i_pathToOutputFile, i_weightsVector, i_suffix):
        l_vectorSize = np.size(i_weightsVector)
               
        l_sourceFile = open(i_pathToOutputFile, 'a')
 
        l_declaration = "  const DATATYPE kernels::optimised::weights"+str(i_suffix)+"["+str(l_vectorSize)+"] __attribute__((aligned(ALIGNMENT))) = { \n"
        # replace all occurrences of 'DATATYPE' with 'float' and 'double', respectively                         
        if(self.m_precision=='SP'):
            l_declaration = re.sub(r'\bDATATYPE\b', 'float', l_declaration)
        elif(self.m_precision=='DP'):
            l_declaration = re.sub(r'\bDATATYPE\b', 'double', l_declaration)
        
        l_content = ""
        for i in range(0, l_vectorSize-1):
            l_content += "    %.16g,\n" % (i_weightsVector[i])
        l_content += "    %.16g\n  }" % (i_weightsVector[-1])
        
        l_sourceFile.write(l_declaration)
        l_sourceFile.write(l_content)
        
        l_sourceFile.close()
        