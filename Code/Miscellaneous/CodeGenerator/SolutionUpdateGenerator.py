#!/bin/env python
##
#
#
#--------------------------------------------------------------
# Generates the code for the element update
# for a specific configuration
#--------------------------------------------------------------
#
#
#
import Backend
from MatmulConfig import MatmulConfig
import FunctionSignatures


class SolutionUpdateGenerator:
    m_config = {}

    # name of generated output file
    m_filename = "solutionUpdate.cpp"

    # total length of vector lduh
    m_vectorLength = -1

    # total length of an soa chunk
    m_chunkSize = -1


    def __init__(self, i_config):
        self.m_config = i_config

        # without padding of lduh, luh
        self.m_chunkSize    = self.m_config['nDof']**self.m_config['nDim']
        self.m_vectorLength = self.m_config['nVar'] * self.m_chunkSize

        # TODO
        # padding of lduh, luh?
        #self.m_chunkSize    = Backend.getSizeWithPadding(self.m_config['nDof']**(self.m_config['nDim']-1))
        #self.m_vectorLength = self.m_config['nVar'] * self.m_chunkSize



    def generateCode(self):
        # write #include's and function signature
        self.__writeHeaderForSolutionUpdate()

        l_file = open(self.m_filename, 'a')

        # we can't use this solution because the toolkit prescribes lduh to be read-only.
        # multiply with inverse of mass matrix (divide by Gaussian weights)
        #l_file.write('  #pragma simd\n')
        #l_file.write('  for(int i=0;i<'+str(self.m_chunkSize)+';i++) {\n')
        #for iVar in range(0, self.m_config['nVar']):
        #    l_offset = iVar * self.m_chunkSize
        #    l_file.write('    lduh['+str(l_offset)+'+i] = 1/kernels::optimised::weights3[i] * lduh['+str(l_offset)+'+i];\n')
        #l_file.write('  }\n\n')
        #
        # sum the contribution to luh
        #l_file.write('  #pragma simd\n')
        #l_file.write('  for(int i=0;i<'+str(self.m_vectorLength)+';i++) {\n')
        #l_file.write('    luh[i] += dt * lduh[i];\n')
        #l_file.write('  }\n\n')

        # We must give the 'restrict' keyword here - otherwise auto-vectorisation fails.
        l_file.write('  #pragma simd\n')
        l_file.write('  for(int i=0;i<'+str(self.m_chunkSize)+';i++) {\n')
        for iVar in range(0, self.m_config['nVar']):
            l_offset = iVar * self.m_chunkSize
            l_file.write('    luh['+str(l_offset)+'+i] += dt/kernels::optimised::weights3[i] * lduh['+str(l_offset)+'+i];\n')
        l_file.write('  }\n\n')

        # write missing closing bracket
        l_file.write('}')
        l_file.close()


    def __writeHeaderForSolutionUpdate(self):
        l_description = '// update the elements \n\n'

        l_includeStatement = '#include "kernels/aderdg/optimised/Kernels.h"\n'   \
                             '#include "kernels/aderdg/optimised/Weights.h"\n\n'

        l_functionSignature = FunctionSignatures.getSolutionUpdateSignature()+" {\n"

        l_sourceFile = open(self.m_filename, 'a')
        l_sourceFile.write(l_description)
        l_sourceFile.write(l_includeStatement)
        l_sourceFile.write(l_functionSignature)
        l_sourceFile.close()


