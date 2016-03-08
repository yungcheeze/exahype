#!/bin/env python

import os
from os.path import join
from os.path import isfile
import subprocess
import errno
from matplotlib.cbook import dedent
from glob import iglob
from shutil import move
import FunctionSignatures
import string
import re


m_architecture           = ''
m_precision              = ''
m_config                 = {}
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
    for l_matmul in i_matmulConfigList:
        l_commandLineArguments =       "dense"  + \
                                 ' ' + m_pathToLibxsmmGenerator+"/"+i_pathToOutputFile + \
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
    
    
def prepareOutputDirectory(i_outputDirectory):
    # create directory for output files if not existing
    try:
        os.makedirs(i_outputDirectory)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
    
    # remove all .cpp and .h files (we are in append mode!) 
    for l_fileName in os.listdir(i_outputDirectory):
        _ , l_ext = os.path.splitext(l_fileName)
        if(l_ext in ['.cpp', '.h']):
            os.remove(i_outputDirectory + "/" + l_fileName)
  


def executeBashCommand(i_command, i_commandLineParameters):
    # usage: executeBashCommand("ls", "-l -a")
    l_bashCommand = i_command + " " + i_commandLineParameters
    l_commandOutput = subprocess.check_output(l_bashCommand, shell=True)
    return l_commandOutput

   
     
def writeIntrinsicsInclude(i_pathToFile):
    l_includeStatement = dedent(  """
                                  #if defined( __SSE3__) || defined(__MIC__) 
                                  #include <immintrin.h>
                                  #endif                                  
                                  """)
    l_includeStatement += "\n\n"
    l_sourceFile = open(i_pathToFile, 'a')
    l_sourceFile.write(l_includeStatement)
    l_sourceFile.close()
     
    
def writeCommonHeader(i_pathToHeaderFile):
    # typically we write to Kernels.h
    l_sourceFile = open(i_pathToHeaderFile, 'a')

    # include guard
    l_sourceFile.write('#ifndef EXAHYPE_KERNELS_OPTIMISED_KERNELS_H_\n'   \
                       '#define EXAHYPE_KERNELS_OPTIMISED_KERNELS_H_\n\n'
                       )

    # global includes
    l_sourceFile.write('#include "tarch/la/Vector.h"\n'                   \
                       '#include "peano/utils/Globals.h"\n\n'
                       )


    # TODO temporary solution
    # necessary till all generic kernels have been replaced with generated ones
    l_sourceFile.write('#define basisSize '+str(m_config['nDof'])+'\n'    \
                       '#define numberOfVariables '+str(m_config['nVar'])+'\n\n')


    if(m_config['nDim']==2):
        pass
        # TODO
    elif(m_config['nDim']==3):
        # global variables/defines
        # TODO: adapt this when padding is introduced
        l_sourceFile.write('#define MbasisSize ' + str(m_config['nDof']) +'\n' \
                           '#define Mvar ' + str(m_config['nVar'])   + '\n'    \
                           '#define Mvar ' + str(m_config['nDim'])   + '\n'    \
                           '#define Mface '+ str(m_config['nDim']*2) + '\n'    \
                           '#define f2p5(var, dim, i, j, k) (var + Mvar*dim + Mvar*Mdim*i + Mvar*Mdim*MbasisSize*j + Mvar*Mdim*MbasisSize*MbasisSize*k)\n' \
                           '#define p2f5(var, dim, i, j, k) (dim*MbasisSize*MbasisSize*MbasisSize*Mvar + Mvar*i + Mvar*MbasisSize*j + Mvar*MbasisSize*MbasisSize*k + var)\n' \
                           '#define f2p4(var, face, a, b) (var + Mvar*face + Mvar*Mface*a + Mvar*Mface*MbasisSize*b)\n' \
                           '#define p2f4(var, face, a, b) (face*MbasisSize*MbasisSize*Mvar + Mvar*a + Mvar*MbasisSize*b + var)\n\n\n'
                           )
    else:
        print("Backend.writeCommonHeader(): nDim not supported")

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
    l_sourceFile.write('#include "kernels/aderdg/optimised/spaceTimePredictor.cpph"\n\n')
    l_sourceFile.write('#include "kernels/aderdg/optimised/riemannSolver.cpph"\n\n')
    
    # close include guard
    l_sourceFile.write('#endif // EXAHYPE_KERNELS_OPTIMISED_KERNELS_H_')
         
    l_sourceFile.close()


def moveGeneratedCppFiles(i_pathToSrc,i_pathToDest):
    l_files = iglob(join(i_pathToSrc, "*.cpp"))
    for l_file in l_files:
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


def setPathToLibxsmmGenerator(i_pathToLibxsmmGenerator):
    global m_pathToLibxsmmGenerator
    m_pathToLibxsmmGenerator = i_pathToLibxsmmGenerator

# -------------------------------------------------------------------
# helpers
# -------------------------------------------------------------------    
def reindentBlock(i_string, i_nSpaces):
    l_string = string.split(i_string, '\n')
    l_string = map(lambda line, ns=i_nSpaces: (' '*ns)+line, l_string)
    l_string = string.join(l_string, '\n')
    return l_string
