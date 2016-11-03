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
# Starting point of the code generator
#
# @note
# requires python3
#
# For a quick test, type
# python Driver.py Euler 5 3 2 nonlinear hsw path/to/libxsmmRepository
#                 NoOtherOption nVar Order 2/3d NoOtherOption Architecture 
#
# for Jenkins this is
# python Driver.py Euler 5 3 2 nonlinear hsw ../../libxsmm --precision=DP
#
# 

import argparse
import CodeGenArgumentParser
from Backend import prepareOutputDirectory
from Backend import moveGeneratedFiles
import Backend
import os
import sys


# --------------------------------------------------------
# Require python3
# --------------------------------------------------------
requiredVersion = (3,0)
currentVersion  = sys.version_info

if(requiredVersion > currentVersion):
    sys.exit('CodeGenerator: Requires Python 3.0 or newer. Abort.')


# --------------------------------------------------------
# Process the command line arguments
# --------------------------------------------------------
l_parser = argparse.ArgumentParser(description='This is the front end of the ExaHyPE code generator.')

l_parser.add_argument('PDE',
                      type=lambda pdeArg: CodeGenArgumentParser.validatePDE(l_parser, pdeArg),
                      help='what example you want to solve')
l_parser.add_argument('numberOfVariables', 
                      type=int, 
                      help='the number of quantities')
l_parser.add_argument('order', 
                      type=int, 
                      help='the order of the approximation polynomial')
l_parser.add_argument('dimension', 
                      type=int, 
                      help='number of dimensions you want to simulate')
l_parser.add_argument('numerics',
                      type=lambda numericsArg: CodeGenArgumentParser.validateNumerics(l_parser, numericsArg),
                      help='linear or nonlinear')
l_parser.add_argument('architecture',
                      type=lambda architectureArg: CodeGenArgumentParser.validateArchitecture(l_parser, architectureArg), 
                      help='the microarchitecture of the target device')
l_parser.add_argument('pathToLibxsmm',
                      type=lambda pathArg: CodeGenArgumentParser.validateLibxsmmGenerator(l_parser, pathArg),
                      help='where to find your local copy of code generator back end "https://github.com/hfp/libxsmm"')
l_parser.add_argument('--precision',
                      type=lambda precisionArg: CodeGenArgumentParser.validatePrecision(l_parser, precisionArg),
                      default='DP',
                      help='SP or DP')


l_commandLineArguments = l_parser.parse_args()

pde                    = l_commandLineArguments.PDE
numberOfVariables      = l_commandLineArguments.numberOfVariables
order                  = l_commandLineArguments.order
dimensions             = l_commandLineArguments.dimension
numerics               = l_commandLineArguments.numerics
architecture           = l_commandLineArguments.architecture
pathToLibxsmmGenerator = l_commandLineArguments.pathToLibxsmm
precision              = l_commandLineArguments.precision

config = { 
           "nVar"              : numberOfVariables,
           "nDof"              : order+1,
           "nDim"              : dimensions
          }

# configure global setup of the code generator
Backend.setArchitecture(architecture)
Backend.setPrecision(precision)
Backend.setConfig(config)
Backend.setNumerics(numerics)
Backend.setPathToLibxsmmGenerator(pathToLibxsmmGenerator)

# clean up output directory
# uncomment when using the Toolkit call chain
#pathToOutputDirectory = "../../Code/ExaHyPE/kernels/aderdg/optimised"
# used for testing as standalone tool
dir = os.path.dirname(__file__)+'/'
pathToOutputDirectory = dir+"../ExaHyPE/kernels/aderdg/optimised"
prepareOutputDirectory(pathToOutputDirectory)

# --------------------------------------------------------
# Now let's generate the compute kernels.
# --------------------------------------------------------

Backend.generateCommonHeader()
Backend.generateComputeKernels()

# --------------------------------------------------------
# Move generated code
# --------------------------------------------------------

# move assembly code
moveGeneratedFiles(pathToLibxsmmGenerator, pathToOutputDirectory)
# move C++ wrapper
moveGeneratedFiles(os.getcwd(), pathToOutputDirectory)

