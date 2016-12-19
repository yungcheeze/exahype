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
# Generates the code for the quadrature weights
# for a specific configuration
#
# @todo TODO
# remove gaussLegendreWeights***

import numpy as np
import Backend
import re

class WeightsGenerator:
    # order of the approximation polynomial
    m_order      = -1

    # number of dimensions we simulate
    m_nDim       = -1

    # quadrature weights mapped onto [0,1]
    m_wGPN       = []

    # name of generated output files
    m_sourceName = "GaussLegendreQuadrature.cpp"
    m_headerName = "GaussLegendreQuadrature.h"

    # spec file-dependent vector with weights combinations
    m_vectors = {}


    def __init__(self, i_config):
        self.m_order     = i_config['nDof']-1
        self.m_nDim      = i_config['nDim']

        # compute the Gauss-Legendre weights
        _, w = np.polynomial.legendre.leggauss(self.m_order+1)
        # map onto [0,1]
        self.m_wGPN = 0.5*w


    def generateCode(self):
        self.__generateHeaderFile()
        self.__generateWeightsCombinations()
        self.__writeToFile()


    def __generateWeightsCombinations(self):
        # We need three kinds of combinations
        # (a) weight = wGPN
        #     suffix :=1
        # (b) aux    = (/ 1.0, wGPN(j), wGPN(k) /)
        #     weight = PRODUCT(aux(1:nDim))
        #     suffix := 2
        # (c) aux    = (/ wGPN(i), wGPN(j), wGPN(k) /)
        #     weight = PRODUCT(aux(1:nDim))
        #     suffix := 3


        # case (1)
        # in any case (2D/3D) we need the ordinary Gauss-Legendre weights
        # in contrast to the generic version, guarantee alignment and pad
        # weight := [wGPN + Pad]
        l_sizeWithoutPadding = np.size(self.m_wGPN)
        l_padWidth           = Backend.getPadWidth(l_sizeWithoutPadding)
        l_weightsVector      = np.pad(self.m_wGPN, (0, l_padWidth), mode='constant')
        self.m_vectors['weights1'] = l_weightsVector


        if(self.m_nDim == 2):
            # case (b)
            # weightsVector is wGPN itself
            # pad weights vector with zeros
            l_sizeWithoutPadding = np.size(self.m_wGPN) 
            l_padWidth           = Backend.getPadWidth(l_sizeWithoutPadding)
            l_weightsVector      = np.pad(self.m_wGPN, (0, l_padWidth), mode='constant')

            self.m_vectors['weights2'] = l_weightsVector

            # case (c)
            # all combinations of two weights, written as an 1D array
            l_weightsVector = np.outer(self.m_wGPN, self.m_wGPN).flatten('F')

            # pad this vector with zeros
            l_sizeWithoutPadding = np.size(l_weightsVector)
            l_padWidth           = Backend.getPadWidth(l_sizeWithoutPadding)
            l_weightsVector      = np.pad(l_weightsVector, (0, l_padWidth), mode='constant')

            self.m_vectors['weights3'] = l_weightsVector

        elif(self.m_nDim == 3):
            # case (b)
            # all combinations of two weights, written as an 1D array
            l_weightsVector = np.outer(self.m_wGPN, self.m_wGPN).flatten('F')

            # pad this vector with zeros
            l_sizeWithoutPadding = np.size(l_weightsVector) 
            l_padWidth           = Backend.getPadWidth(l_sizeWithoutPadding)
            l_weightsVector      = np.pad(l_weightsVector, (0, l_padWidth), mode='constant')

            self.m_vectors['weights2'] = l_weightsVector

            # case (c)
            # all combination of three weights, written as an 1D array
            l_weightsVector = np.kron(np.outer(self.m_wGPN, self.m_wGPN), self.m_wGPN).flatten('F')

            # pad this vector with zeros
            l_sizeWithoutPadding = np.size(l_weightsVector)
            l_padWidth           = Backend.getPadWidth(l_sizeWithoutPadding)
            l_weightsVector      = np.pad(l_weightsVector, (0, l_padWidth), mode='constant')

            self.m_vectors['weights3'] = l_weightsVector

        else:
            print("WeightsGenerator.__generateWeightsCombinations(): nDim not supported")



    def __generateHeaderFile(self):
        l_sourceFile = open(self.m_headerName, 'a')
        # mostly copied from Dominic's GaussLegendreQuadrature.h
        l_includeGuard = '#ifndef _EXAHYPE_KERNELS_ADERDG_OPTIMISED_WEIGHTS_H_\n'   \
                         '#define _EXAHYPE_KERNELS_ADERDG_OPTIMISED_WEIGHTS_H_\n\n'
        l_sourceFile.write(l_includeGuard)
        l_sourceFile.write('#include <set>\n\n')
        l_sourceFile.write('namespace kernels { \n'                                                  \
                           'namespace aderdg {\n'                                                    \
                           'namespace optimised {\n\n'                                               \
                           'void initGaussLegendreNodesAndWeights(const std::set<int>& orders);\n'   \
                           'void freeGaussLegendreNodesAndWeights(const std::set<int>& orders);\n\n' \
                           'extern double **gaussLegendreNodes;\n'                                   \
                           'extern double **gaussLegendreWeights;\n'                                 \
                           'extern double *weights1;\n'                                              \
                           'extern double *weights2;\n'                                              \
                           'extern double *weights3;\n}\n}\n}\n')
        # close include guard
        l_sourceFile.write('#endif /* _EXAHYPE_KERNELS_ADERDG_OPTIMISED_WEIGHTS_H_ */')
        l_sourceFile.close()


    def __writeToFile(self):
        l_sourceFile = open(self.m_sourceName, 'a')
        l_sourceFile.write('#include "kernels/aderdg/optimised/'+self.m_headerName+'"\n' \
                           '#include <mm_malloc.h> //g++\n\n')
        l_sourceFile.write('double** kernels::aderdg::optimised::gaussLegendreWeights;\n'  \
                           'double** kernels::aderdg::optimised::gaussLegendreNodes;\n'    \
                           'double* kernels::aderdg::optimised::weights1;\n'               \
                           'double* kernels::aderdg::optimised::weights2;\n'               \
                           'double* kernels::aderdg::optimised::weights3;\n\n')

        l_sourceFile.write('void kernels::aderdg::optimised::freeGaussLegendreNodesAndWeights(const std::set<int>& orders) {\n')
        for weightsVector in self.m_vectors:
            l_sourceFile.write('  _mm_free('+str(weightsVector)+');\n')
        l_sourceFile.write('  constexpr int MAX_ORDER=9;\n\n'              \
                           '  for (int i = 0; i < MAX_ORDER + 1; i++) {\n' \
                           '    delete [] gaussLegendreNodes[i];\n'        \
                           '    delete [] gaussLegendreWeights[i];\n'      \
                           '  }\n\n'                                       \
                           '  delete [] gaussLegendreNodes;\n'             \
                           '  delete [] gaussLegendreWeights;\n')
        l_sourceFile.write('}\n\n')

        l_sourceFile.write('void kernels::aderdg::optimised::initGaussLegendreNodesAndWeights(const std::set<int>& orders) {\n')
        for weightsVector in self.m_vectors:
            l_elemCount = self.m_vectors[weightsVector].size
            l_sourceFile.write('  '+str(weightsVector)+' = (double *) _mm_malloc(sizeof(double)*'+ \
                               str(l_elemCount)+', ALIGNMENT);\n')

        for weightsVector in self.m_vectors:
            for idx in range(0, self.m_vectors[weightsVector].size):
                l_sourceFile.write('  '+str(weightsVector) + '['+str(idx)+'] = %.12e' % self.m_vectors[weightsVector].item(idx)+';\n')

        # for compatibility with generic code
        l_sourceFile.write('  constexpr int MAX_ORDER=9;\n'\
                            '\n'\
                            '  gaussLegendreNodes = new double* [MAX_ORDER + 1];\n'\
                            '  gaussLegendreWeights = new double* [MAX_ORDER + 1];\n'\
                            '\n'\
                            '  for (int i = 0; i < MAX_ORDER + 1; i++) {\n'\
                            '    gaussLegendreNodes[i] = new double[i + 1];\n'\
                            '    gaussLegendreWeights[i] = new double[i + 1];\n'\
                            '  }\n'\
                            '\n'\
                            '  gaussLegendreWeights[0][0] = 1.0000000000000000;\n'\
                            '  gaussLegendreNodes[0][0] = 0.5000000000000000;\n'\
                            '\n'\
                            '  gaussLegendreWeights[1][0] = 0.5000000000000000;\n'\
                            '  gaussLegendreWeights[1][1] = 0.5000000000000000;\n'\
                            '  gaussLegendreNodes[1][0] = 0.2113248654051871;\n'\
                            '  gaussLegendreNodes[1][1] = 0.7886751345948129;\n'\
                            '\n'\
                            '  gaussLegendreWeights[2][0] = 0.2777777777777778;\n'\
                            '  gaussLegendreWeights[2][1] = 0.4444444444444444;\n'\
                            '  gaussLegendreWeights[2][2] = 0.2777777777777778;\n'\
                            '  gaussLegendreNodes[2][0] = 0.1127016653792583;\n'\
                            '  gaussLegendreNodes[2][1] = 0.5000000000000000;\n'\
                            '  gaussLegendreNodes[2][2] = 0.8872983346207417;\n'\
                            '\n'\
                            '  gaussLegendreWeights[3][0] = 0.1739274225687273;\n'\
                            '  gaussLegendreWeights[3][1] = 0.3260725774312732;\n'\
                            '  gaussLegendreWeights[3][2] = 0.3260725774312732;\n'\
                            '  gaussLegendreWeights[3][3] = 0.1739274225687273;\n'\
                            '  gaussLegendreNodes[3][0] = 0.0694318442029737;\n'\
                            '  gaussLegendreNodes[3][1] = 0.3300094782075719;\n'\
                            '  gaussLegendreNodes[3][2] = 0.6699905217924281;\n'\
                            '  gaussLegendreNodes[3][3] = 0.9305681557970262;\n'\
                            '\n'\
                            '  gaussLegendreWeights[4][0] = 0.1184634425280948;\n'\
                            '  gaussLegendreWeights[4][1] = 0.239314335249683;\n'\
                            '  gaussLegendreWeights[4][2] = 0.2844444444444443;\n'\
                            '  gaussLegendreWeights[4][3] = 0.239314335249683;\n'\
                            '  gaussLegendreWeights[4][4] = 0.1184634425280948;\n'\
                            '  gaussLegendreNodes[4][0] = 0.04691007703066802;\n'\
                            '  gaussLegendreNodes[4][1] = 0.2307653449471584;\n'\
                            '  gaussLegendreNodes[4][2] = 0.5000000000000000;\n'\
                            '  gaussLegendreNodes[4][3] = 0.7692346550528415;\n'\
                            '  gaussLegendreNodes[4][4] = 0.9530899229693319;\n'\
                            '\n'\
                            '  gaussLegendreWeights[5][0] = 0.0856622461895845;\n'\
                            '  gaussLegendreWeights[5][1] = 0.1803807865240695;\n'\
                            '  gaussLegendreWeights[5][2] = 0.2339569672863459;\n'\
                            '  gaussLegendreWeights[5][3] = 0.2339569672863459;\n'\
                            '  gaussLegendreWeights[5][4] = 0.1803807865240695;\n'\
                            '  gaussLegendreWeights[5][5] = 0.0856622461895845;\n'\
                            '  gaussLegendreNodes[5][0] = 0.03376524289842397;\n'\
                            '  gaussLegendreNodes[5][1] = 0.1693953067668678;\n'\
                            '  gaussLegendreNodes[5][2] = 0.3806904069584015;\n'\
                            '  gaussLegendreNodes[5][3] = 0.6193095930415985;\n'\
                            '  gaussLegendreNodes[5][4] = 0.8306046932331322;\n'\
                            '  gaussLegendreNodes[5][5] = 0.966234757101576;\n'\
                            '\n'\
                            '  gaussLegendreWeights[6][0] = 0.06474248308443538;\n'\
                            '  gaussLegendreWeights[6][1] = 0.1398526957446382;\n'\
                            '  gaussLegendreWeights[6][2] = 0.1909150252525592;\n'\
                            '  gaussLegendreWeights[6][3] = 0.2089795918367344;\n'\
                            '  gaussLegendreWeights[6][4] = 0.1909150252525592;\n'\
                            '  gaussLegendreWeights[6][5] = 0.1398526957446382;\n'\
                            '  gaussLegendreWeights[6][6] = 0.06474248308443538;\n'\
                            '  gaussLegendreNodes[6][0] = 0.02544604382862076;\n'\
                            '  gaussLegendreNodes[6][1] = 0.1292344072003028;\n'\
                            '  gaussLegendreNodes[6][2] = 0.2970774243113014;\n'\
                            '  gaussLegendreNodes[6][3] = 0.5000000000000000;\n'\
                            '  gaussLegendreNodes[6][4] = 0.7029225756886985;\n'\
                            '  gaussLegendreNodes[6][5] = 0.8707655927996972;\n'\
                            '  gaussLegendreNodes[6][6] = 0.9745539561713792;\n'\
                            '\n'\
                            '  gaussLegendreWeights[7][0] = 0.05061426814518821;\n'\
                            '  gaussLegendreWeights[7][1] = 0.1111905172266871;\n'\
                            '  gaussLegendreWeights[7][2] = 0.1568533229389437;\n'\
                            '  gaussLegendreWeights[7][3] = 0.1813418916891809;\n'\
                            '  gaussLegendreWeights[7][4] = 0.1813418916891809;\n'\
                            '  gaussLegendreWeights[7][5] = 0.1568533229389437;\n'\
                            '  gaussLegendreWeights[7][6] = 0.1111905172266871;\n'\
                            '  gaussLegendreWeights[7][7] = 0.05061426814518821;\n'\
                            '  gaussLegendreNodes[7][0] = 0.01985507175123186;\n'\
                            '  gaussLegendreNodes[7][1] = 0.1016667612931866;\n'\
                            '  gaussLegendreNodes[7][2] = 0.2372337950418355;\n'\
                            '  gaussLegendreNodes[7][3] = 0.4082826787521751;\n'\
                            '  gaussLegendreNodes[7][4] = 0.5917173212478249;\n'\
                            '  gaussLegendreNodes[7][5] = 0.7627662049581645;\n'\
                            '  gaussLegendreNodes[7][6] = 0.8983332387068134;\n'\
                            '  gaussLegendreNodes[7][7] = 0.9801449282487682;\n'\
                            '\n'\
                            '  gaussLegendreWeights[8][0] = 0.04063719418078751;\n'\
                            '  gaussLegendreWeights[8][1] = 0.09032408034742861;\n'\
                            '  gaussLegendreWeights[8][2] = 0.1303053482014677;\n'\
                            '  gaussLegendreWeights[8][3] = 0.1561735385200013;\n'\
                            '  gaussLegendreWeights[8][4] = 0.1651196775006297;\n'\
                            '  gaussLegendreWeights[8][5] = 0.1561735385200013;\n'\
                            '  gaussLegendreWeights[8][6] = 0.1303053482014677;\n'\
                            '  gaussLegendreWeights[8][7] = 0.09032408034742861;\n'\
                            '  gaussLegendreWeights[8][8] = 0.04063719418078751;\n'\
                            '  gaussLegendreNodes[8][0] = 0.01591988024618696;\n'\
                            '  gaussLegendreNodes[8][1] = 0.08198444633668212;\n'\
                            '  gaussLegendreNodes[8][2] = 0.1933142836497048;\n'\
                            '  gaussLegendreNodes[8][3] = 0.3378732882980955;\n'\
                            '  gaussLegendreNodes[8][4] = 0.5000000000000000;\n'\
                            '  gaussLegendreNodes[8][5] = 0.6621267117019045;\n'\
                            '  gaussLegendreNodes[8][6] = 0.8066857163502952;\n'\
                            '  gaussLegendreNodes[8][7] = 0.9180155536633179;\n'\
                            '  gaussLegendreNodes[8][8] = 0.984080119753813;\n'\
                            '\n'\
                            '  gaussLegendreWeights[9][0] = 0.03333567215434358;\n'\
                            '  gaussLegendreWeights[9][1] = 0.07472567457529024;\n'\
                            '  gaussLegendreWeights[9][2] = 0.1095431812579912;\n'\
                            '  gaussLegendreWeights[9][3] = 0.1346333596549983;\n'\
                            '  gaussLegendreWeights[9][4] = 0.1477621123573766;\n'\
                            '  gaussLegendreWeights[9][5] = 0.1477621123573766;\n'\
                            '  gaussLegendreWeights[9][6] = 0.1346333596549983;\n'\
                            '  gaussLegendreWeights[9][7] = 0.1095431812579912;\n'\
                            '  gaussLegendreWeights[9][8] = 0.07472567457529024;\n'\
                            '  gaussLegendreWeights[9][9] = 0.03333567215434358;\n'\
                            '  gaussLegendreNodes[9][0] = 0.01304673574141413;\n'\
                            '  gaussLegendreNodes[9][1] = 0.06746831665550773;\n'\
                            '  gaussLegendreNodes[9][2] = 0.1602952158504878;\n'\
                            '  gaussLegendreNodes[9][3] = 0.2833023029353764;\n'\
                            '  gaussLegendreNodes[9][4] = 0.4255628305091844;\n'\
                            '  gaussLegendreNodes[9][5] = 0.5744371694908156;\n'\
                            '  gaussLegendreNodes[9][6] = 0.7166976970646236;\n'\
                            '  gaussLegendreNodes[9][7] = 0.8397047841495122;\n'\
                            '  gaussLegendreNodes[9][8] = 0.9325316833444923;\n'\
                            '  gaussLegendreNodes[9][9] = 0.9869532642585859;\n')

        l_sourceFile.write("}\n\n")
        l_sourceFile.close()












