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


m_architecture           = ''
m_precision              = ''
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
    
    # remove all .cpp files (we are in append mode!) 
    for fileName in os.listdir(i_outputDirectory):
        if fileName.endswith(".cpp"):
            os.remove(i_outputDirectory + "/" + fileName)    


def executeBashCommand(i_command, i_commandLineParameters):
    # usage: executeBashCommand("ls", "-l -a")
    l_bashCommand = i_command + " " + i_commandLineParameters
    l_commandOutput = subprocess.check_output(l_bashCommand, shell=True)
    return l_commandOutput


def validateLibxsmmGenerator(i_pathToLibxsmm):
    l_pathToLibxsmmGenerator = i_pathToLibxsmm + "/bin/libxsmm_gemm_generator"
    return isfile(l_pathToLibxsmmGenerator)
     
     
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
     
    
def writeCommonKernelsHeaderInclude(i_pathToFile):
    l_includeStatement = dedent(  """
                                  #include "kernels/aderdg/optimised/Kernels.h"                                 
                                  """)
    l_includeStatement += "\n\n"
    l_sourceFile = open(i_pathToFile, 'a')
    l_sourceFile.write(l_includeStatement)
    l_sourceFile.close()


def moveGeneratedCppFiles(i_pathToSrc,i_pathToDest):
    l_files = iglob(join(i_pathToSrc, "*.cpp"))
    for l_file in l_files:
        if(isfile(l_file)):
            move(l_file, i_pathToDest)


def setArchitecture(i_architecture):
    global m_architecture 
    m_architecture = i_architecture


def setPrecision(i_precision):
    global m_precision
    m_precision = i_precision
    FunctionSignatures.setPrecision(i_precision)


def setPathToLibxsmmGenerator(i_pathToLibxsmmGenerator):
    global m_pathToLibxsmmGenerator
    m_pathToLibxsmmGenerator = i_pathToLibxsmmGenerator