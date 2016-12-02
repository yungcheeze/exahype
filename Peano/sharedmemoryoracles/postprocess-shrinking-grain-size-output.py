import sys
import re
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


def processMeasurement(adapter):
  searchPattern = "adapter=" + str(adapter) + ","

  searchPattern = "adapter-number=" + str(adapter) 
  line          = inputFile.readline()
  while not line=="" and not re.search( searchPattern, line ):
    line = inputFile.readline()
  
  htmlOverview.write( "<table border=\"1\">" );
  htmlOverview.write( "<tr><td><b>Method</b></td><td><b>Grain Size</b></td><td><b>Maximum problem size</b></td><td><b>Accuracy</b></td><td><b>Remarks</b></td></tr>" );

  line = inputFile.readline()
  while not re.search( "end OracleForOnePhaseWithShrinkingGrainSize", line):
    methodTrace        = line.split( "=")[0]
    biggestProblemSize = line.split( "=")[1].split( "," )[0]
    grainSize          = line.split( "=")[1].split( "," )[1]
    previousGrainSize  = line.split( "=")[1].split( "," )[2]
    accuracy           = line.split( "eps=")[1].split( "," )[0]
    
    isNotScaling       = int(grainSize)==0
    isStillSearching   = float(accuracy)>1e-8
    
    htmlOverview.write( "<tr>" );
    htmlOverview.write( "<td>" + methodTrace + "</td>" );
    if isNotScaling and isStillSearching:
      htmlOverview.write( "<td bgcolor=\"yellow\">" + grainSize + "</td>" );
    elif isNotScaling:
      htmlOverview.write( "<td bgcolor=\"red\">" + grainSize + "</td>" );
    else:
      htmlOverview.write( "<td bgcolor=\"green\">" + grainSize + "</td>" );
    htmlOverview.write( "<td>" + biggestProblemSize + "</td>" );
    if isStillSearching:
      htmlOverview.write( "<td>" + accuracy + "</td>" );
    else:
      htmlOverview.write( "<td bgcolor=\"red\">" + accuracy + "</td>" );
    
    htmlOverview.write( "<td>" );
    if isNotScaling:
      htmlOverview.write( "Seems not to scale. " );
    if isStillSearching:
      htmlOverview.write( "Still searching for accurate data. " );
    htmlOverview.write( "</td>" );
    htmlOverview.write( "</tr>" );
    line = inputFile.readline()  

  htmlOverview.write( "</table>" );
      
  return;

#
#
#   main
# ========
#
#    
if (len(sys.argv)!=2):
  print "usage: python ../postprocess-shrinking-grain-size-output.py outputfile"
  quit()

htmlOverview = open( sys.argv[1] + ".html",  "w" )
htmlOverview.write( "<h1>" + sys.argv[1] + "</h1>" );

inputFile = open(sys.argv[1], "r" )

line                  = inputFile.readline()
totalNumberOfOracles  = int( line.split("=")[1] )
line                  = inputFile.readline()
oraclesForSteering    = int( line.split("=")[1] )

htmlOverview.write( "Number of oraces=" + str(totalNumberOfOracles) + " (incl. oracles required for repository/algorithm steering)" );
htmlOverview.write( "<br />" );
htmlOverview.write( "Oracles required for repository steering=" + str(oraclesForSteering) );

htmlOverview.write( "<h3>Table of content:</h3>" );
htmlOverview.write( "<ul>" );
for adapter in range(oraclesForSteering,totalNumberOfOracles):
  htmlOverview.write( "<li><a href=\"#adapter-" + str(adapter) + "\">Adapter " + str(adapter) + "</a></li>" );
htmlOverview.write( "</ul>" );

htmlOverview.write( "<p>Empty adapter sections imply that this adapter is not used by the code. The adapter order requals the adapter order in the specification file. Only those program phases that are actually used are also displayed.</p>" );

for adapter in range(oraclesForSteering,totalNumberOfOracles):
  htmlOverview.write( "<h3 id=\"adapter-" + str(adapter) + "\">Adapter " + str(adapter) + "</h3>" );
  processMeasurement(adapter)
