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

    # total length of vectors luh, lduh
    m_vectorLength = -1

    # total length of an soa chunk
    m_chunkSize = -1

    # number of spatials dofs
    m_nDof = -1

    def __init__(self, i_config):
        self.m_config = i_config

        # total number of spatial dofs
        self.m_nDof = self.m_config['nDof']**self.m_config['nDim']

        # length luh, lduh
        self.m_vectorLength = self.m_config['nVar']*self.m_nDof

        # version: soa format
        # without padding of lduh, luh
        #self.m_chunkSize    = self.m_config['nDof']**self.m_config['nDim']
        #self.m_vectorLength = self.m_config['nVar'] * self.m_chunkSize
        #
        # alternative schemes for soa format. padding of lduh, luh?
        #self.m_chunkSize    = Backend.getSizeWithPadding(self.m_config['nDof']**(self.m_config['nDim']-1))
        #self.m_vectorLength = self.m_config['nVar'] * self.m_chunkSize
        #
        # alternatively
        #self.m_chunkSize    = self.m_config['nDof']**self.m_config['nDim']
        #self.m_vectorLength = Backend.getSizeWithPadding(self.m_chunkSize)



    def generateCode(self):
        # write #include's and function signature
        self.__writeHeaderForSolutionUpdate()

        l_sourceFile = open(self.m_filename, 'a')

        # version: soa format
        # We must give the 'restrict' keyword here - otherwise auto-vectorisation fails.
        #l_file.write('  #pragma simd\n')
        #l_file.write('  for(int i=0;i<'+str(self.m_chunkSize)+';i++) {\n')
        #for iVar in range(0, self.m_config['nVar']):
        #    l_offset = iVar * self.m_chunkSize
        #    l_file.write('    luh['+str(l_offset)+'+i] += dt/kernels::optimised::weights3[i] * lduh['+str(l_offset)+'+i];\n')
        #l_file.write('  }\n\n')

        # version: aos format (use this for compatibility with export functions, etc. for the time being)
        # need temporary memory because the toolkit prescribes lduh to be read-only.
        l_sourceFile.write('  double tmp['+str(self.m_vectorLength)+'] __attribute__((aligned(ALIGNMENT)));\n\n')
        # multiply with inverse of mass matrix (divide by Gaussian weights)
        # this part is not vectorised by the compiler. It's anyway only a temporary solution.
        l_sourceFile.write('#pragma simd\n'\
                           '  for(int i=0;i<'+str(self.m_nDof)+';i++) {\n')
        for iVar in range(0, self.m_config['nVar']):
            l_sourceFile.write('    tmp[i*'+str(self.m_config['nVar'])+'+'+str(iVar)+'] = lduh[i*'+str(self.m_config['nVar'])+'+'+str(iVar)+']/kernels::weights3[i];\n')
        l_sourceFile.write('  }\n\n')

        # sum the contribution to lduh
        # this part is vectorised by the compiler
        l_sourceFile.write('#pragma simd\n'\
                           '  for(int i=0;i<'+str(self.m_vectorLength)+';i++) {\n'\
                           '    luh[i] += dt * tmp[i];\n'\
                           '  }\n')

        # write missing closing bracket
        l_sourceFile.write('}')
        l_sourceFile.close()


    def __writeHeaderForSolutionUpdate(self):
        l_description = '// update the elements \n\n'

        l_includeStatement = '#include "kernels/aderdg/optimised/Kernels.h"\n'   \
                             '#include "kernels/aderdg/optimised/GaussLegendreQuadrature.h"\n\n'

        l_functionSignature = FunctionSignatures.getSolutionUpdateSignature()+" {\n"

        l_sourceFile = open(self.m_filename, 'a')
        l_sourceFile.write(l_description)
        l_sourceFile.write(l_includeStatement)
        l_sourceFile.write(l_functionSignature)
        l_sourceFile.close()


