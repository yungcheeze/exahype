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


import TemplatingUtils


class KernelsHeaderGenerator:
    m_context = {}
    m_args = {'useV2' : False}

    # name of generated output file
    m_filename = 'Kernels.h'

    
    def __init__(self, i_context, i_args):
        self.m_context = i_context
        self.m_args = i_args


    def generateCode(self):
        self.m_context['solverNamespace'] = self.m_context['solverName'].split('::')[0]
        self.m_context['solverClass'] = self.m_context['solverName'].split('::')[1]
        if(self.m_args['useV2']): 
            TemplatingUtils.renderAsFile('KernelsV2_h.template', self.m_filename, self.m_context)
        else:
            TemplatingUtils.renderAsFile('Kernels_h.template', self.m_filename, self.m_context)

