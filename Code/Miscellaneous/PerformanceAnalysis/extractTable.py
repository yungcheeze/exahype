import argparse
from argparse import RawTextHelpFormatter

import os
import runtimeParser
import hpclib 


def createTable(path,prefix,postfix,MaxNodes):
  print "write table for prefix " + prefix + " and postfix " + postfix + ". Max count searched for is " + str(MaxNodes)
  outFilename   = prefix + postfix + ".table"
  outFile       = open( outFilename, "w" )
  adapters      = []
  for ranks in range(1,MaxNodes):
    inputFileName = path + "/" + prefix + str(ranks) + postfix
    print "- search for " + inputFileName + " ... ",
    if os.path.isfile(inputFileName):
      print "found "

      times    = runtimeParser.parse_adapter_times(inputFileName)
      maxLevel = runtimeParser.max_level(inputFileName)
      newLine  = "" # We have to build up the line completely but only write it 
                    # in the end. Otherwise we might get incomplete lines 
      
      if len(adapters)==0 and len(times)>0:
        outFile.write( "ranks/threads/nodes " )
        outFile.write( "& " )
        outFile.write( "max level " )
        for t in times:
          outFile.write( " & " )
          outFile.write( t )
          adapters.append( t )
        outFile.write( "& Total " )
        outFile.write( "\n" )
        
      if len(adapters)>0:
        newLine += str(ranks)
        newLine += " & "  
        newLine += str(maxLevel) 
        totalIterationCount = 0
        totalRuntime        = 0
        for adapter in adapters:
          newLine += " & " 
          newLine += str(times[adapter]['n']) 
          totalIterationCount = totalIterationCount + times[adapter]['n']
          newLine += " & " 
          newLine += str(times[adapter]['usertime']) 
          totalRuntime = totalRuntime + times[adapter]['usertime']
      
        newLine +=  " & " 
        newLine += str(totalIterationCount) 
        newLine += " & " 
        newLine += str(totalRuntime) 
        
      if newLine!="":
        outFile.write( newLine + "\n" ) 
    else:
      print "not available "

########################################################################
# START OF THE PROGRAM
########################################################################
# For the meaning of the command line arguments, see the documentation
# of the function 'plot_multithreading_adapter_scaling' above.
help = '''
Creates a speedup plot based on Peano output files with specific file naming pattern.
If multiple adapters are specified, then the cumulative user times are used to compute the speedup.
\n\n
Sample usage:\n\n
python plot-cluster-speedup.py -path tmp
\n\n
The code searches for files with the pattern prefixNpostfix, where N is the number of compute nodes.
'''

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('-path',required=True,help="Directory containing the Peano output files.")
parser.add_argument('-prefix',required=True,help="Prefix of the files to be parsed. You can add \"\" if no prefix is used.")
parser.add_argument('-postfix',required=True,help="Postfix of the files to be parsed. You can add \"\" if no postfix is used. If you postfix starts with a minus, please embed it into \\\' quotes.")
parser.add_argument('-maxnodes',required=True,help="Maximum number of nodes.")

args = parser.parse_args();

createTable(args.path,args.prefix,args.postfix,int(args.maxnodes))

