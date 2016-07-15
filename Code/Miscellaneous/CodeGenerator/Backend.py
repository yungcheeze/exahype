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
# This file is pivotal to the code generator. It manages 
# the internal decisions about padding, triggers the code
# generation of the solver kernels and calls the back end 
# assembly code generation.
#

import os
from os.path import join
from os.path import isfile
import subprocess
import errno
from glob import iglob
from shutil import move
import FunctionSignatures
import SpaceTimePredictorGenerator
import VolumeIntegralGenerator
import SurfaceIntegralGenerator
import RiemannGenerator
import SolutionUpdateGenerator
import StableTimeStepSizeGenerator
import WeightsGenerator
import DGMatrixGenerator
import string
import re


m_architecture           = ''
m_precision              = ''
m_config                 = {}
m_numerics               = ''
m_pathToLibxsmmGenerator = ''
m_simdWidth              =  {'SP':  {'noarch' : 1,
                                     'wsm'    : 4,
                                     'snb'    : 8,
                                     'hsw'    : 8,
                                     'knc'    : 16,
                                     'knl'    : 16 },
                             'DP': { 'noarch' : 1,
                                     'wsm'    : 2,
                                     'snb'    : 4,
                                     'hsw'    : 4,
                                     'knc'    : 8,
                                     'knl'    : 8 }
                            }


def executeLibxsmmGenerator(i_commandLineParameters):
    l_bashCommand = m_pathToLibxsmmGenerator + '/libxsmm_gemm_generator ' + i_commandLineParameters
    subprocess.call(l_bashCommand.split())


def generateAssemblerCode(i_pathToOutputFile,
                          i_matmulConfigList):
    l_pathToAsmFile = os.path.splitext(i_pathToOutputFile)[0]+'.c'
    for l_matmul in i_matmulConfigList:
        l_commandLineArguments =       "dense"  + \
                                 ' ' + m_pathToLibxsmmGenerator+"/"+l_pathToAsmFile + \
                                 ' ' + l_matmul.baseroutinename + \
                                 ' ' + str(l_matmul.M) + \
                                 ' ' + str(l_matmul.N) + \
                                 ' ' + str(l_matmul.K) + \
                                 ' ' + str(l_matmul.LDA) + \
                                 ' ' + str(l_matmul.LDB) + \
                                 ' ' + str(l_matmul.LDC) + \
                                 ' ' + str(l_matmul.alpha) + \
                                 ' ' + str(l_matmul.beta) + \
                                 ' ' + str(l_matmul.alignment_A) + \
                                 ' ' + str(l_matmul.alignment_C) + \
                                 ' ' + m_architecture + \
                                 ' ' + l_matmul.prefetchStrategy+ \
                                 ' ' + m_precision 
        executeLibxsmmGenerator(l_commandLineArguments)


def getSizeWithPadding(i_sizeWithoutPadding: int) -> int:
    l_simdSize        = m_simdWidth[m_precision][m_architecture]
    l_sizeWithPadding = l_simdSize * int((i_sizeWithoutPadding+(l_simdSize-1))/l_simdSize)
    return l_sizeWithPadding


def getPadWidth(i_sizeWithoutPadding: int) -> int:
    l_simdSize        = m_simdWidth[m_precision][m_architecture]
    l_sizeWithPadding = l_simdSize * int((i_sizeWithoutPadding+(l_simdSize-1))/l_simdSize)
    l_padWidth        = l_sizeWithPadding - i_sizeWithoutPadding
    return l_padWidth


def prepareOutputDirectory(i_outputDirectory):
    # create directory for output files if not existing
    try:
        os.makedirs(i_outputDirectory)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    # remove all .cpp, .cpph and .h files (we are in append mode!)
    for l_fileName in os.listdir(i_outputDirectory):
        _ , l_ext = os.path.splitext(l_fileName)
        if(l_ext in ['.cpp', '.cpph', '.c', '.h']):
            os.remove(i_outputDirectory + "/" + l_fileName)



def executeBashCommand(i_command, i_commandLineParameters):
    # usage: executeBashCommand("ls", "-l -a")
    l_bashCommand = i_command + " " + i_commandLineParameters
    l_commandOutput = subprocess.check_output(l_bashCommand, shell=True)
    return l_commandOutput



def generateCommonHeader():
    # name of generated output file
    l_filename = "Kernels.h"
    l_sourceFile = open(l_filename, 'a')

    # include guard
    l_sourceFile.write('#ifndef EXAHYPE_KERNELS_OPTIMISED_KERNELS_H_\n'   \
                       '#define EXAHYPE_KERNELS_OPTIMISED_KERNELS_H_\n\n'
                       )

    # global includes
    l_sourceFile.write('#include "tarch/la/Vector.h"\n'                   \
                       '#include "peano/utils/Globals.h"\n\n'
                       )

    # nested namespaces
    l_sourceFile.write('namespace kernels {\n'        )
    l_sourceFile.write('  namespace aderdg {\n'       )
    l_sourceFile.write('    namespace optimised {\n\n')

    # we are inside of a nested namespace and add extra spaces to our function signatures
    l_indentation = 6

    # fetch function signatures
    l_functionList = []
    l_functionList.append(FunctionSignatures.getPicardLoopSignature(m_config['nDim']))
    l_functionList.append(FunctionSignatures.getPredictorSignature())
    l_functionList.append(FunctionSignatures.getExtrapolatorSignature())
    l_functionList.append(FunctionSignatures.getSolutionUpdateSignature())
    l_functionList.append(FunctionSignatures.getVolumeIntegralSignature())
    l_functionList.append(FunctionSignatures.getSurfaceIntegralSignature())
    l_functionList.append(FunctionSignatures.getInitialConditionSignature())
    l_functionList.append(FunctionSignatures.getSolutionAdjustmentSignature())
    l_functionList.append(FunctionSignatures.getRiemannSolverSignature())
    l_functionList.append(FunctionSignatures.getStableTimeStepSizeSignature())

    # declare c++ functions in header file
    for l_functionSignature in l_functionList:
        # we are already inside the namespace, so we cut off the namespace substring
        l_functionSignature = re.sub('kernels::aderdg::optimised::','',l_functionSignature)
        # fix indentation
        l_functionSignature = reindentBlock(l_functionSignature, l_indentation)
        # write function declarations one after another
        l_sourceFile.write(l_functionSignature+";\n\n")


    # closing brackets of namespace
    l_sourceFile.write('    }\n')
    l_sourceFile.write('  }\n'  )
    l_sourceFile.write('}\n\n'  )

    # include template functions
    l_sourceFile.write('#include "kernels/aderdg/optimised/solutionAdjustment.cpph"\n\n')
    l_sourceFile.write('#include "kernels/aderdg/optimised/stableTimeStepSize.cpph"\n\n')
    if(m_numerics == 'nonlinear'):
        l_sourceFile.write('#include "kernels/aderdg/optimised/picard.cpph"\n\n')
    else:
        l_sourceFile.write('#include "kernels/aderdg/optimised/cauchyKovalewski.cpph"\n\n')
    l_sourceFile.write('#include "kernels/aderdg/optimised/riemannSolver.cpph"\n\n')

    # close include guard
    l_sourceFile.write('#endif // EXAHYPE_KERNELS_OPTIMISED_KERNELS_H_')

    l_sourceFile.close()


def generateComputeKernels():
    spaceTimePredictorGenerator = SpaceTimePredictorGenerator.SpaceTimePredictorGenerator(m_config, m_numerics)
    spaceTimePredictorGenerator.generateCode()
    volumeIntegralGenerator = VolumeIntegralGenerator.VolumeIntegralGenerator(m_config, m_numerics)
    volumeIntegralGenerator.generateCode()
    surfaceIntegralGenerator = SurfaceIntegralGenerator.SurfaceIntegralGenerator(m_config, m_numerics)
    surfaceIntegralGenerator.generateCode()
    riemannGenerator = RiemannGenerator.RiemannGenerator(m_config, m_numerics, m_precision)
    riemannGenerator.generateCode()
    solutionUpdateGenerator = SolutionUpdateGenerator.SolutionUpdateGenerator(m_config)
    solutionUpdateGenerator.generateCode()
    stableTimeStepSizeGenerator = StableTimeStepSizeGenerator.StableTimeStepSizeGenerator(m_config)
    stableTimeStepSizeGenerator.generateCode()
    weightsGenerator = WeightsGenerator.WeightsGenerator(m_config)
    weightsGenerator.generateCode()
    dgMatrixGenerator = DGMatrixGenerator.DGMatrixGenerator(m_config, m_numerics)
    dgMatrixGenerator.generateCode()



def moveGeneratedFiles(i_pathToSrc,i_pathToDest):
    l_fileTypes = ('*.h', '*.cpp', '*.c', '*.cpph')
    l_fileList = []
    for l_file in l_fileTypes:
        l_fileList.extend(iglob(i_pathToSrc+"/"+l_file))

    for l_file in l_fileList:
        if(isfile(l_file)):
            move(l_file, i_pathToDest)

# -------------------------------------------------------------------
# Class variables 
# -------------------------------------------------------------------
def setArchitecture(i_architecture):
    global m_architecture 
    m_architecture = i_architecture


def setPrecision(i_precision):
    global m_precision
    m_precision = i_precision
    FunctionSignatures.setPrecision(i_precision)


def setConfig(i_config):
    global m_config
    m_config = i_config

def setNumerics(i_numerics):
    global m_numerics
    m_numerics = i_numerics

def setPathToLibxsmmGenerator(i_pathToLibxsmmGenerator):
    global m_pathToLibxsmmGenerator
    m_pathToLibxsmmGenerator = i_pathToLibxsmmGenerator

# -------------------------------------------------------------------
# helpers
# -------------------------------------------------------------------
def reindentLine(i_line, i_nSpaces=8):
    return (' ' * i_nSpaces) + i_line

def reindentBlock(i_string, i_nSpaces):
    l_stringList = i_string.split('\n')
    l_stringList = list(map(lambda line, ns=i_nSpaces: reindentLine(line, ns), l_stringList))
    l_string = "\n".join(l_stringList)
    return l_string