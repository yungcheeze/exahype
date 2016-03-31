import sys
import re
import pylab 
import os


Symbol = [ 
  "s", "o", ">", 
  "<", "^", "v" 
]
Colour = [
  "#ff0000", "#00ff00", "#0000ff",
  "#ffff00", "#ff00ff", "#00ffff"
]
AlphaValue=0.5


def processMeasurement(adapter,phase):
  try:
    searchPattern = "adapter=" + str(adapter) + ", phase=" + str(phase) + ","
    substring    = line.split( searchPattern )[1]
    substring    = substring.split( "adapter" )[0] # everything left of the remaining string
  except:
    return ""

  #print "found data for adapter=" + str(adapter) + ", phase=" + str(phase) + ": " + substring
  
  name                = substring.split("name=")[1].split(",")[0] 
  studiedProblemSizes = substring.split("no-of-problem-sizes-studied=")[1].split(",")[0] 

  htmlOverview.write( "name=" + name + "<br />" );
  htmlOverview.write( "adapter no=" + str(adapter) + ", phase no=" + str(phase) + "<br />" );
  htmlOverview.write( "no of studied problem sizes=" + studiedProblemSizes + "<br />" );

  pylab.clf()
  currentColourAndSymbol = 0
  for problemSize in range(0,int(studiedProblemSizes)):
    problemSizeSubString = substring.split("problem-size=")[problemSize+1]
    problemSize          = problemSizeSubString.split(",")[0] 
    grainSizes           = problemSizeSubString.split("no-of-grain-sizes-studied=")[1].split(":")[0] 

    xValues = []
    yValues = []
    zValues = []
    
    for measurement in problemSizeSubString.split( "(grain-size=" )[1:]:
      grainSize  = float(measurement.split( "," )[0])
      time       = float(measurement.split( "(" )[1].split( "," )[0])
      deviation  = float(measurement.split( "deviation=" )[1].split( ")" )[0])
      xValues.append( grainSize )
      yValues.append( time )
      zValues.append( deviation  )

    pylab.plot(xValues, yValues, "-" + Symbol[currentColourAndSymbol],  markersize=10, color=Colour[currentColourAndSymbol], label=str(problemSize) )
    pylab.plot(xValues, zValues,       Symbol[currentColourAndSymbol],  markersize=12, color=Colour[currentColourAndSymbol], alpha=AlphaValue )
    currentColourAndSymbol = currentColourAndSymbol+1
    
  outputFileName = "adapter-" + str(adapter) + "-phase-" + str(phase)
  print "write file " + outputFileName
    
  symbolAndColourCounter = 0
  
  try:     
    pylab.grid(True)
    pylab.ylabel('[t]=s')
    pylab.yscale( 'log' )
    if (xValues[-1]>100):
      pylab.xscale( 'symlog' )
    pylab.xlabel('grain size')
    pylab.legend(fontsize=9, framealpha=0.5)
    pylab.legend(loc='upper left',framealpha=0.5)
  
    pylab.savefig( outputFileName + ".png")
    pylab.savefig( outputFileName + ".pdf")
    htmlOverview.write( "<img src=\"" + outputFileName + ".png\" /> <br />" );
  except:
    print "Unexpected error:", sys.exc_info()[0]
    return ""

  return name





#
#
#   main
# ========
#
#    
if (len(sys.argv)!=2):
  print "usage: python ../postprocess-sampling-output.py outputfile"
  #print "  no-of-adapter  number of adapters you have in your project (see your specification file)"
  #print "  no-of-phases   number of code phases that are tuned via an oracle. 19 by default (cmp OracleForOnePhase)"
  quit()

htmlOverview = open( sys.argv[1] + ".html",  "w" )
htmlOverview.write( "<h1>" + sys.argv[1] + "</h1>" );

inputFile = open(sys.argv[1], "r" )

line                  = inputFile.readline()
totalNumberOfOracles  = int( line.split("=")[1] )
line                  = inputFile.readline()
numberOfAdapters      = int( line.split("=")[1] )
line                  = inputFile.readline()
adaptersForSteering   = int( line.split("=")[1] )
line                  = inputFile.readline()
noOfMethodsCalling    = int( line.split("=")[1] )

htmlOverview.write( "Total number of oracles=" + str(totalNumberOfOracles) );
htmlOverview.write( "<br />" );
htmlOverview.write( "Number of adapters=" + str(numberOfAdapters) + " (incl. adapters required for algorithm steering)" );
htmlOverview.write( "<br />" );
htmlOverview.write( "Adapters required for repository steering=" + str(adaptersForSteering) );
htmlOverview.write( "<br />" );
htmlOverview.write( "Number of methods calling=" + str(noOfMethodsCalling) + " (max)");

htmlOverview.write( 
  "<p>The solid lines present timings for particular problem sizes plus grain sizes. The semitransparent symbols denote the standard deviation belonging to the measurements. Each pair of symbols plus symbols with lines studies one particular problem size over multiple grain sizes. The legend gives the problem size.</p>" 
);

htmlOverview.write( "<h3>Table of content:</h3>" );
htmlOverview.write( "<ul>" );
for adapter in range(0,numberOfAdapters-adaptersForSteering):
  htmlOverview.write( "<li><a href=\"#adapter-" + str(adapter) + "\">Adapter " + str(adapter) + "</a></li>" );
htmlOverview.write( "</ul>" );

line                  = inputFile.readline()   # The remainder is all written into one line
for adapter in range(0,numberOfAdapters-adaptersForSteering):
  htmlOverview.write( "<h3 id=\"adapter-" + str(adapter) + "\">Adapter " + str(adapter) + "</h3>" );
  for phase in range(0,noOfMethodsCalling):
    processMeasurement(adapter,phase)
