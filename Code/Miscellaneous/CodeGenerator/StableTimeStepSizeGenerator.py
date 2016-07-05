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
# Generates the code for the computation of the time step
# It's simplistic and just picks the 2D or 3D version
# of the generic kernels.
#

import FunctionSignatures
import Backend

class StableTimeStepSizeGenerator:
    m_config = {}

    # order of the approximation polynomial
    m_order = -1

    # number of degrees of freedom
    m_nDof = -1

    # 2D, 3D
    m_nDim  = -1

    # name of generated output file
    m_filename = "stableTimeStepSize.cpph"

    # maximum CFL numbers, copied from generic code base
    m_PNPM = [1.0,   0.33,  0.17, 0.1,  0.069, 0.045, 0.038, 0.03, 0.02, 0.015]

    # CFL number corrected by PNPM condition
    m_CFL = -1

    # number of quantities
    m_nVar = -1


    def __init__(self, i_config):
        self.m_order  = i_config['nDof']-1
        self.m_nDim   = i_config['nDim']
        self.m_CFL    = 0.9 * self.m_PNPM[self.m_order]
        self.m_nVar   = i_config['nVar']
        self.m_nDof   = i_config['nDof']
        self.m_config = i_config


    def generateCode(self):
        self.__writeHeader()

        l_sourceFile = open(self.m_filename, 'a')

        # temporary arrays for the eigenvalues and assembly of quantities
        l_paddedVectorLength = Backend.getSizeWithPadding(self.m_nVar)

        #-----------------------------------------------------------------------------------
        # soa version. later on we probably want to revive this.
        #-----------------------------------------------------------------------------------
        ## without padding
        #l_chunkSize = self.m_nDof ** self.m_nDim
        #l_vectorLength = self.m_nVar * l_chunkSize

        ## with padding
        ##l_chunkSize = Backend.getSizeWithPadding(self.m_nDof ** self.m_nDim)
        ##l_vectorLength = self.m_nVar * l_chunkSize

        ## number of real Dofs per chunk without any padded entries
        #l_varsPerChunk = self.m_nDof**self.m_nDim

        #l_sourceFile.write('  double lambda['+str(l_paddedVectorLength)+'] __attribute__((aligned(ALIGNMENT)));\n')
        #l_sourceFile.write('  double contiguousVars['+str(l_paddedVectorLength)+'] __attribute__((aligned(ALIGNMENT)));\n\n')
        #l_sourceFile.write('  double dt = std::numeric_limits<double>::max();\n')

        #l_sourceFile.write('  for(int i=0;i<'+str(l_varsPerChunk)+';i++) {\n')
        #for iVar in range(0, self.m_nVar):
            #l_sourceFile.write('    contiguousVars['+str(iVar)+'] = luh[i+'+str(iVar*l_chunkSize)+'];\n')

        #l_sourceFile.write('\n')
        #l_sourceFile.write('    double denominator = 0.0;\n')
        #l_sourceFile.write('    for(int d=0;d<'+str(self.m_nDim)+';d++) {\n')
        #l_sourceFile.write('      PDEEigenvalues(&contiguousVars[0], d, &lambda[0]);\n\n')
        #l_sourceFile.write('      double maxEigenvalue = 0.0;\n')
        ## process without(!) padding
        #l_sourceFile.write('      for (int ivar = 0; ivar < '+str(self.m_nVar)+'; ivar++) {\n')
        #l_sourceFile.write('        maxEigenvalue = std::max(fabs(lambda[ivar]), maxEigenvalue);\n')
        #l_sourceFile.write('      }\n')
        #l_sourceFile.write('      denominator += maxEigenvalue / dx[d];\n')
        #l_sourceFile.write('    }\n\n')
        #l_sourceFile.write('    dt = std::min(dt, ' + str(self.m_CFL) + '/denominator);\n')
        #l_sourceFile.write('  }\n')
        #l_sourceFile.write('  return dt;\n')
        #l_sourceFile.write('}')
        #l_sourceFile.close()

        #-----------------------------------------------------------------------------------
        # aos version. temporary solution. be compatible with generic code base
        # ugly and slow. directly copied from the generic code base. 
        #-----------------------------------------------------------------------------------
        l_sourceFile.write('  double lambda['+str(l_paddedVectorLength)+'] __attribute__((aligned(ALIGNMENT)));\n')
        l_sourceFile.write('  double dt = std::numeric_limits<double>::max();\n')
        if(self.m_config['nDim']==2):
            l_sourceFile.write('  for(int ii=0;ii<'+str(self.m_config['nDof'])+';ii++) {\n'\
                               '    for(int jj=0;jj<'+str(self.m_config['nDof'])+';jj++) {\n'\
                               '      const int nodeIndex = ii + '+str(self.m_config['nDof'])+' * jj;\n'\
                               '      const int dofStartIndex = nodeIndex * '+str(self.m_config['nVar'])+';\n'\
                               '      double denominator = 0.0;\n'\
                               '      for (int d = 0; d < 2; d++) {\n'\
                               '        PDEEigenvalues(&luh[dofStartIndex], d, lambda);\n\n'\
                               '        double maxEigenvalue = 0.0;\n'\
                               '        for (int ivar = 0; ivar < '+str(self.m_config['nVar'])+'; ivar++) {\n'\
                               '          maxEigenvalue = std::max(fabs(lambda[ivar]), maxEigenvalue);\n'\
                               '        }\n'\
                               '        denominator += maxEigenvalue / dx[d];\n'\
                               '      }\n\n'\
                               '      dt = std::min(dt, ' + str(self.m_CFL) + '/denominator);\n'\
                               '    }\n'\
                               '  }\n')
        elif(self.m_config['nDim']==3):
            l_sourceFile.write('  for(int ii=0;ii<'+str(self.m_config['nDof'])+';ii++) {\n'\
                               '    for(int jj=0;jj<'+str(self.m_config['nDof'])+';jj++) {\n'\
                               '      for(int kk=0;kk<'+str(self.m_config['nDof'])+';kk++) {\n'\
                               '        const int nodeIndex = ii + '+str(self.m_config['nDof'])+' * jj + '+str(self.m_config['nDof']**2)+' * kk;\n'\
                               '        const int dofStartIndex = nodeIndex * '+str(self.m_config['nVar'])+';\n'\
                               '        double denominator = 0.0;\n'\
                               '        for (int d = 0; d < 3; d++) {\n'\
                               '          PDEEigenvalues(&luh[dofStartIndex], d, lambda);\n\n'\
                               '          double maxEigenvalue = 0.0;\n'\
                               '          for (int ivar = 0; ivar < '+str(self.m_config['nVar'])+'; ivar++) {\n'\
                               '            maxEigenvalue = std::max(fabs(lambda[ivar]), maxEigenvalue);\n'\
                               '          }\n'\
                               '          denominator += maxEigenvalue / dx[d];\n'\
                               '        }\n\n'\
                               '        dt = std::min(dt, ' + str(self.m_CFL) + '/denominator);\n'\
                               '      }\n'\
                               '    }\n'\
                               '  }\n')

        l_sourceFile.write('  return dt;\n')
        l_sourceFile.write('}')
        l_sourceFile.close()


    def __writeHeader(self):
        l_includeStatement = '#include "kernels/aderdg/optimised/Kernels.h"\n' \
                             '#include <limits>\n\n'

        l_functionSignature = FunctionSignatures.getStableTimeStepSizeSignature()+" {\n"

        l_sourceFile = open(self.m_filename, 'a')
        l_sourceFile.write(l_includeStatement)
        l_sourceFile.write(l_functionSignature)
        l_sourceFile.close()