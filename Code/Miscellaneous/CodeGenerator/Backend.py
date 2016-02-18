#!/bin/env python

import os
import subprocess
from os import remove
import errno

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
            remove(i_outputDirectory + "/" + fileName)    



def executeBashCommand(i_command, i_commandLineParameters):
    # usage: executeBashCommand("ls", "-l -a")
    l_bashCommand = i_command + " " + i_commandLineParameters
    l_commandOutput = subprocess.check_output(l_bashCommand, shell=True)
    return l_commandOutput


def validateLibxsmmGenerator(i_pathToLibxsmm):
    l_pathToLibxsmmGenerator = i_pathToLibxsmm + "/bin/libxsmm_gemm_generator"
    return os.path.isfile(l_pathToLibxsmmGenerator)
     
    