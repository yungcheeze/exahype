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
# Generates the code for the mapping of the DG polynomial
# onto the [0,1] and calls the user-defined function.
# It's simplistic and just picks the 2D or 3D version
# of the generic kernels.
#

import os
from jinja2 import Template

import Backend


class SolutionAdjustmentGenerator:
    m_context = {}

    # name of generated output file
    m_filename = "solutionAdjustment.cpph"


    def __init__(self, i_context):
        self.m_context = i_context


    def generateCode(self):
        dir = os.path.dirname(__file__)+'/'
        self.m_context['order'] = self.m_context['nDof']-1
        self.m_context['nDimPad'] = Backend.getSizeWithPadding(self.m_context['nDim'])
        
        with open(dir+'templates/solutionAdjustment.template', 'r') as tmp:
            template = Template(tmp.read())
            with open(self.m_filename, 'w') as out:
                out.write(template.render(self.m_context))
