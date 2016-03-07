#!/bin/env python
##
#
#
#--------------------------------------------------------------
# Starting point of the code generator
#--------------------------------------------------------------
#
# For a quick test, type
# python Driver.py MyEulerSolver 5 3 2 hsw path/to/libxsmmRepository
#
# 


import argparse
import CodeGenArgumentParser
from SpaceTimePredictorGenerator import SpaceTimePredictorGenerator
from Backend import prepareOutputDirectory
from Backend import moveGeneratedCppFiles
import Backend
import os
       

# --------------------------------------------------------
# Process the command line arguments
# --------------------------------------------------------
l_parser = argparse.ArgumentParser(description='This is the frontend of the ExaHyPE code generator.')

l_parser.add_argument('solverName', 
                      type=str,
                      help='the namespace')
l_parser.add_argument('numberOfVariables', 
                      type=int, 
                      help='the number of quantities')
l_parser.add_argument('order', 
                      type=int, 
                      help='the order of the approximation polynomial')
l_parser.add_argument('dimension', 
                      type=int, 
                      help='number of dimensions you want to simulate')
l_parser.add_argument('architecture',
                      type=lambda architectureArg: CodeGenArgumentParser.validateArchitecture(l_parser, architectureArg), 
                      help='the microarchitecture of the target device')
l_parser.add_argument('pathToLibxsmm',
                      type=lambda pathArg: CodeGenArgumentParser.validateLibxsmmGenerator(l_parser, pathArg),
                      help='where to find your local copy of code generator backend "https://github.com/hfp/libxsmm"')
l_parser.add_argument('--precision',
                      type=lambda precisionArg: CodeGenArgumentParser.validatePrecision(l_parser, precisionArg),
                      default='DP',
                      help='SP or DP')


l_commandLineArguments = l_parser.parse_args()

solverName             = l_commandLineArguments.solverName
numberOfVariables      = l_commandLineArguments.numberOfVariables
order                  = l_commandLineArguments.order
dimensions             = l_commandLineArguments.dimension
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
Backend.setPathToLibxsmmGenerator(pathToLibxsmmGenerator)

# clean up output directory
pathToOutputDirectory = "../../ExaHyPE/kernels/aderdg/optimised"
prepareOutputDirectory(pathToOutputDirectory)

# --------------------------------------------------------
# Now let's generate the compute kernels.
# --------------------------------------------------------

Backend.writeCommonHeader("Kernels.h")
# TODO move Kernels.h to directory kernels/aderdg/optimised 

spaceTimePredictorGenerator = SpaceTimePredictorGenerator(config)
spaceTimePredictorGenerator.generateCode()

# move assembler code
moveGeneratedCppFiles(pathToLibxsmmGenerator, pathToOutputDirectory)
# move C++ wrapper
moveGeneratedCppFiles(os.getcwd(), pathToOutputDirectory)


