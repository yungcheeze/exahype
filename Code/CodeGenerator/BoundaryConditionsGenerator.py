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


import os
from jinja2 import Template

import Backend


class BoundaryConditionsGenerator:
    m_context = {}

    # linear/nonlinear
    m_type   = ""

    # name of generated output file
    m_filenameRoot = "boundaryConditions"

    
    
    def __init__(self, i_config, i_numerics):
        self.m_context = i_config
        self.m_type   = i_numerics


    def generateCode(self):
        dir = os.path.dirname(__file__)
        self.m_context['blockSize'] = self.m_context['nDofPad'] if self.m_context['nDim'] == 2 else Backend.getSizeWithPadding(self.m_context['nDof'] * self.m_context['nDof'])
        self.m_context['iVar_range_0_nVar'] = range(0, self.m_context['nVar'])

        with open(os.path.join(dir,'templates/boundaryConditions.template'), 'r') as tmp:
            template = Template(tmp.read(), trim_blocks=True)
            with open(self.m_filenameRoot+'.cpph', 'w') as out:
                out.write(template.render(self.m_context))
                    
