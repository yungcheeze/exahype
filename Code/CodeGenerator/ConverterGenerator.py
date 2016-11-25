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

class ConverterGenerator:
    m_context = {}

    # linear/nonlinear
    m_type   = ""

    # name of generated output file
    m_filenameRoot = "converter"



    def __init__(self, i_config, i_numerics):
        self.m_context = i_config
        self.m_type   = i_numerics
        


    def generateCode(self):
        dir = os.path.dirname(__file__)
        self.m_context['bndBlockSize'] = self.m_context['nDofPad'] if self.m_context['nDim'] == 2 else getSizeWithPadding(self.m_context['nDof'] * self.m_context['nDof'])
    
        with open(os.path.join(dir,'templates/converter_h.template'), 'r') as tmp:
            template = Template(tmp.read(), trim_blocks=True)
            with open(self.m_filenameRoot+'.h', 'w') as out:
                out.write(template.render(self.m_context))
                
        with open(os.path.join(dir,'templates/converter_cpp.template'), 'r') as tmp:
            template = Template(tmp.read(), trim_blocks=True)
            with open(self.m_filenameRoot+'.cpp', 'w') as out:
                out.write(template.render(self.m_context))



















