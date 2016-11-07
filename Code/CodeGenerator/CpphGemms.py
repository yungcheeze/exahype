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
from MatmulConfig import MatmulConfig

class CpphGemms:
    m_context = {}

    # linear/nonlinear
    m_type   = ""

    # name of generated output file
    m_filenameRoot = "CpphGemms"



    def __init__(self, i_config, i_numerics):
        self.m_context = i_config
        self.m_type   = i_numerics


    def generateCode(self):
        dir = os.path.dirname(__file__)+'/'
        self.m_context['gemm_a_b_c']  = 'gemm_'+str(self.m_context['nVar'])+'_'+str(self.m_context['nDof'])+'_'+str(self.m_context['nDof'])
    
        if(self.m_type == 'linear'):
            pass
        else:
            with open(dir+'templates/CpphGemms_h.template', 'r') as tmp:
                template = Template(tmp.read())
                with open(self.m_filenameRoot+'.h', 'w') as out:
                    out.write(template.render(self.m_context))
                    
            with open(dir+'templates/CpphGemms_cpp.template', 'r') as tmp:
                template = Template(tmp.read())
                with open(self.m_filenameRoot+'.cpp', 'w') as out:
                    out.write(template.render(self.m_context))



















