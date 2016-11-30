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
# Generates the code for the surface integral
# for a specific configuration
#


import Backend
import TemplatingUtils


class SurfaceIntegralGenerator:
    m_context = {}

    # linear/nonlinear
    m_type   = ""

    # name of generated output file
    m_filename = "surfaceIntegral.cpp"


    def __init__(self, i_context, i_numerics):
        self.m_context = i_context
        self.m_type    = i_numerics


    def generateCode(self):
        self.m_context['bndBlockSize'] = self.m_context['nDofPad'] if self.m_context['nDim'] == 2 else Backend.getSizeWithPadding(self.m_context['nDof'] * self.m_context['nDof'])
        self.m_context['bndFaceSize'] = self.m_context['nVar'] * self.m_context['bndBlockSize']
        self.m_context['nDofPowDimMinOne'] = self.m_context['nDof'] if self.m_context['nDim'] == 2 else (self.m_context['nDof'] * self.m_context['nDof'])
    
        if(self.m_type == 'linear'):
           pass
        else:
            TemplatingUtils.renderAsFile('surfaceIntegralNonLinear_cpp.template', self.m_filename, self.m_context)






