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


class SpaceTimePredictorGenerator:
    m_context = {}

    # name of generated output file
    m_filename_picard       = 'picard.cpph'
    m_filename_predictor    = 'predictor.cpp'
    m_filename_extrapolator = 'extrapolatedPredictor.cpp'

    
    def __init__(self, i_config):
        self.m_context = i_config


    def generateCode(self):
        self.m_context['bndBlockSize'] = self.m_context['nDofPad'] if self.m_context['nDim'] == 2 else Backend.getSizeWithPadding(self.m_context['nDof'] * self.m_context['nDof'])
        self.m_context['nDof3D'] = 1 if self.m_context['nDim'] == 2 else self.m_context['nDof']
        
        TemplatingUtils.renderAsFile('spaceTimePredictor_picard_cpph.template', self.m_filename_picard, self.m_context)
        TemplatingUtils.renderAsFile('spaceTimePredictor_predictor_cpp.template', self.m_filename_predictor, self.m_context)
        TemplatingUtils.renderAsFile('spaceTimePredictor_extrapolator_cpp.template', self.m_filename_extrapolator, self.m_context)

