#! /usr/bin/env python

# This file is part of the ExaHyPE project.
# Copyright (c) 2016  http://exahype.eu
# All rights reserved.
# 
# The project has received funding from the European Union's Horizon 
# 2020 research and innovation programme under grant agreement
# No 671698. For copyrights and licensing, please consult the webpage.
# 
# Released under the BSD 3 Open Source License.
# For the full license text, see LICENSE.txt

import os
from solutionAdjustmentTest import prepareSolutionAdjustmentTest

regenerateOptimisedKernel = True
dimension = 3
nvar = 5
order = 3
arch = 'hsw'
alignment = 16


prepareSolutionAdjustmentTest(dimension, nvar, order, arch, alignment, regenerateOptimisedKernel)

#Remove __pycache__
print("Removing __pycache__ with: rm -r __pycache__")
os.system('rm -r __pycache__')