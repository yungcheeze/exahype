#!/bin/env python
##
#
#
#--------------------------------------------------------------
# Starting point of the code generator
#--------------------------------------------------------------
#
# For a quick test, type
# python Driver.py MyEulerSolver 5 3 2 nonlinear hsw path/to/libxsmmRepository
#
# for Jenkins this is
# python Driver.py MyEulerSolver 5 3 2 nonlinear hsw ../../../libxsmm --precision=DP
#
# 


import argparse
import CodeGenArgumentParser
from SpaceTimePredictorGenerator import SpaceTimePredictorGenerator
from RiemannGenerator import RiemannGenerator
from Backend import prepareOutputDirectory
from Backend import moveGeneratedFiles
import Backend
import os
#from WeightsGenerator import WeightsGenerator


# --------------------------------------------------------
# Process the command line arguments
# --------------------------------------------------------
l_parser = argparse.ArgumentParser(description='This is the front end of the ExaHyPE code generator.')

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

solverName             = l_commandLineArguments.solverName
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
Backend.setPathToLibxsmmGenerator(pathToLibxsmmGenerator)

# clean up output directory
# uncomment when using the Toolkit call chain
#pathToOutputDirectory = "../../Code/ExaHyPE/kernels/aderdg/optimised"
# used for testing as standalone tool
pathToOutputDirectory = "../../ExaHyPE/kernels/aderdg/optimised"
prepareOutputDirectory(pathToOutputDirectory)

# --------------------------------------------------------
# Now let's generate the compute kernels.
# --------------------------------------------------------

Backend.generateCommonHeader()
spaceTimePredictorGenerator = SpaceTimePredictorGenerator(config)
spaceTimePredictorGenerator.generateCode()
riemannGenerator = RiemannGenerator(config, numerics, precision)
riemannGenerator.generateCode()

# for testing
#weightsGenerator = WeightsGenerator(config, precision)
#weightsGenerator.generateCode()
# TODO move Weights.h to output directory
# extend the moveGeneratedCppFiles() to also work for .h files


# move assembler code
moveGeneratedFiles(pathToLibxsmmGenerator, pathToOutputDirectory)
# move C++ wrapper
moveGeneratedFiles(os.getcwd(), pathToOutputDirectory)

