#!/bin/env python

import os
from os.path import join
from os.path import isfile
import subprocess
import errno
from matplotlib.cbook import dedent
from distlib._backport.shutil import move
from glob import iglob



def executeLibxsmmGenerator(i_pathToLibxsmmGenerator,
                            i_commandLineParameters):
    l_bashCommand = i_pathToLibxsmmGenerator + '/libxsmm_gemm_generator ' + i_commandLineParameters
    subprocess.call(l_bashCommand.split())
    
    
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
