#!/bin/env python
##
#
#
#--------------------------------------------------------------
# Starting point of the code generator
#--------------------------------------------------------------
#
# For a quick test, type
# python Driver.py MyEulerSolver 5 3 2 hsw
#
# 


import argparse
from Backend import validateLibxsmmGenerator
from SpaceTimePredictorGenerator import SpaceTimePredictorGenerator
from Backend import prepareOutputDirectory


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
                      type=str, 
                      help='the microarchitecture of the target device')
l_parser.add_argument('--precision',
                      type=str,
                      default='DP',
                      help='SP or DP')
l_parser.add_argument('--pathToLibxsmm', 
                      help='where to find your local copy of code generator backend "https://github.com/hfp/libxsmm"')

l_commandLineArguments = l_parser.parse_args()

solverName        = l_commandLineArguments.solverName
numberOfVariables = l_commandLineArguments.numberOfVariables
order             = l_commandLineArguments.order
dimensions        = l_commandLineArguments.dimension
architecture      = l_commandLineArguments.architecture


#
# process optional arguments
#
if(vars(l_commandLineArguments)['precision']=='DP'):
    precision    = 'DP'
elif(vars(l_commandLineArguments)['precision']=='SP'):
    precision    = 'SP'
else:
    print("Unknown precision specified. Continue with double precision")
    precision    = 'DP'

# TODO make pathToLibxsmm argument obligatory
if(l_commandLineArguments.pathToLibxsmm == None):
    # just for testing :-)
    l_pathToLibxsmm        = "/home/schwara/Documents/workspace/ga63cad_workspace/libxsmm/"
    pathToLibxsmmGenerator = "/home/schwara/Documents/workspace/ga63cad_workspace/libxsmm/bin"
else:
    l_pathToLibxsmm        = l_commandLineArguments.pathToLibxsmm
    pathToLibxsmmGenerator = l_pathToLibxsmm + "/bin"
    
      
if(not validateLibxsmmGenerator(l_pathToLibxsmm)):
    print("Can't find the code generator of libxsmm. Did you run 'make generator' to build the library?")
    exit()
    
       
config = { 
           "solverName"        : solverName,
           "nVar"              : numberOfVariables,
           "order"             : order,
           "nDim"              : dimensions,
           "architecture"      : architecture,
           "precision"         : precision
          }



# clean up output directory
pathToOutputDirectory = "../../ExaHyPE/kernels/aderdg/optimised"
prepareOutputDirectory(pathToOutputDirectory)

# TODO verify that Kernels.h exists
# create from scratch?

# --------------------------------------------------------
# Now let's generate the compute kernels.
# --------------------------------------------------------

spaceTimePredictorGenerator = SpaceTimePredictorGenerator(config)
spaceTimePredictorGenerator.generateCode()
