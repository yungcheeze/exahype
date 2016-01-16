from PIL import Image, ImageDraw, ImageFont

import csv
import sys
import re
import subprocess



def getRuntime(fileName,adapterNumber,trace,problemSize,grainSize):
  result = 0
  try:
    searchStringHeader =  str(adapterNumber) + "th adapter, method " + trace + ", problem size " + problemSize
    file = open(fileName, "r" )
    data = "   "
    while ( not data=="" and not re.search( searchStringHeader, data ) ):
      data = file.readline()
    searchStringGrainSize = "grain-size: " + str(grainSize)
    data = file.readline()
    while ( 
      not data=="" and 
      not re.search( searchStringGrainSize, data ) 
      and
      not re.search( "runtimes for", data )
    ):
      data = file.readline()
    if (not data=="" and re.search( searchStringGrainSize, data ) ):
      result = data.split(",")[0].split("(")[-1]
  finally:
    return result


def getStudiedProblemSizes(fileName,adapterNumber,methodTrace):
  result = []
  try:
    searchString =  str(adapterNumber) + "th adapter, method " + methodTrace + ", problem size"
    file = open(fileName, "r" )
    data = "   "
    while ( not data=="" ):
      if ( re.search( searchString, data ) ):
        appendEntry = data.split(searchString)[1].split(")")[0].strip()
        result.append(appendEntry)
      data = file.readline()
  finally:
    return result


def prepareGnuplotFile(scriptDir,replace_dict):
  lines = open( scriptDir + '/postprocess-grain-size-sampling.gnuplot.template').read() 
 
  for search, replace in replace_dict.items():
    lines = lines.replace(search, replace)
 
  outFile = open('postprocess-grain-size-sampling.gnuplot', 'w')
  outFile.write(lines) 
  outFile.close()




traces = [
  "load-vertices",
  "load-vertices-on-regular-stationary-grid",
  "load-vertices-on-irregular-stationary-grid",
  "load-vertices-on-stationary-with-parallel-boundary",
  "store-vertices",
  "store-vertices-on-irregular-stationary-grid",
  "store-vertices-on-regular-stationary-grid",
  "set-counter",
  "enter-cell-and-load-sub-cells",
  "leave-cell-and-store-sub-cells",
  "call-enter-cell-and-initialise-enumerators-on-regular-stationary-grid",
  "call-touch-first-time-on-regular-stationary-grid",
  "call-touch-last-time-on-regular-stationary-grid",
  "call-enter-cell-on-regular-stationary-grid",
  "call-leave-cell-on-regular-stationary-grid",
  "pipeline-ascend-task",
  "pipeline-descend-task",
  "split-load-vertices-task-on-regular-stationary-grid",
  "split-store-vertices-task-on-regular-stationary-grid",
  "ascend-on-regular-stationary-grid",
  "descend-on-regular-stationary-grid"
]



if len(sys.argv) != 4:
    print 'usage: python script.py <filename> <adapter-number> <gnuplot script dir>'
    sys.exit(0)

fileName      = sys.argv[1]
adapterNumber = sys.argv[2]
scriptDir     = sys.argv[3]



 
 
for trace in traces:
  outFilenamePrefix = "adapter-" + str(adapterNumber) + "_" + trace
  outFilename = outFilenamePrefix + ".table"
  outFile     = open(outFilename, "w" )
  print( "write file " + outFilename );
  
  valueMax = 0
  
  problemSizes = getStudiedProblemSizes(fileName,adapterNumber,trace)
  outFile.write( "grain-size\\problem-size" )
  for problemSize in problemSizes:
    outFile.write( " | " )
    outFile.write( problemSize )
  outFile.write( "\n" )

  if len(problemSizes)>0:
    print 'will analyse ' + problemSizes[-1] + ' possible different grain sizes'
    for grainSize in range(0,int(problemSizes[-1])):
      foundDatum = False
      for problemSizeCounter,problemSize in enumerate(problemSizes):
        t = getRuntime(fileName,adapterNumber,trace,problemSize,grainSize)
        if (not t==0):
          if (not foundDatum):
            foundDatum = True
            outFile.write( str(grainSize) )
            for i in range(0,problemSizeCounter): 
              outFile.write( " | ?0 " )
          outFile.write( " | " )
          
          value = float(t)/float(problemSize)
          if valueMax<value:
             valueMax = value
          outFile.write( str( value ) )
        if t==0 and foundDatum:
          outFile.write( " | ?0 " )
      if (foundDatum):
        outFile.write( "\n" )
  
    outFile.close() 
  
    replace_dict = {
       'TITLE':    trace,   
       'DATAFILE': outFilename,
       'OUTFILE':  outFilenamePrefix,
       'MAXVALUE': str(valueMax),
       'MAXPROBLEMSIZE': problemSizes[-1]
    }
    for problemSizeCounter,problemSize in enumerate(problemSizes):
      replace_dict[ 'PROBLEMSIZE' + str(problemSizeCounter) ] = str(problemSize)
    
    replace_dict[ 'DEVICE']    = "postscript eps enhanced color"
    replace_dict[ 'EXTENSION'] = "eps"
    prepareGnuplotFile(scriptDir,replace_dict)  
    subprocess.call('gnuplot "postprocess-grain-size-sampling.gnuplot" 2>/dev/null', shell=True);

    replace_dict[ 'DEVICE']    = "pdf enhanced color"
    replace_dict[ 'EXTENSION'] = "pdf"
    prepareGnuplotFile(scriptDir,replace_dict)  
    subprocess.call('gnuplot "postprocess-grain-size-sampling.gnuplot" 2>/dev/null', shell=True);

    replace_dict[ 'DEVICE']    = "png size 1200,900"
    replace_dict[ 'EXTENSION'] = "png"
    prepareGnuplotFile(scriptDir,replace_dict)  
    subprocess.call('gnuplot "postprocess-grain-size-sampling.gnuplot" 2>/dev/null', shell=True);
  else:
    print 'no results available for algorithm phase ' + trace