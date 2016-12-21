import sys
import re
import os
from enum import Enum


Symbol = [ 
  "s", "o", ">", 
  "<", "^", "v" 
]
Colour = [
  "#ff0000", "#00ff00", "#0000ff",
  "#ffff00", "#ff00ff", "#00ffff"
]


class Analysis(Enum):
  Nop                = 0,
  NoSerialRuntimeYet = 1,
  SeemsNotToScale    = 2,
  MightScale       = 3,
  DoesNotScale       = 4,
  Scales             = 5


def processMeasurement(adapter):
  searchPattern = "adapter=" + str(adapter) + ","

  searchPattern = "adapter-number=" + str(adapter) 
  line          = inputFile.readline()
  runtime       = 0
  while not line=="" and not re.search( searchPattern, line ):
    if re.search( "total-runtime=", line ):
      runtime = line.split("=")[1]
    line = inputFile.readline()
  
  htmlOverview.write( "<p>Total runtime=" + runtime + "</p>" );
  htmlOverview.write( "<table border=\"1\">" );
  htmlOverview.write( "<tr><td><b>Method</b></td><td><b>Maximum problem size</b></td><td><b>Grain size</b></td><td><b>Search delta</b></td><td><b>Accuracy</b></td><td><b>Total serial runtime</b></td><td><b>Serial runtime</b></td><td><b>Current runtime</b></td><td><b>Est. speedup</b></td><td><b>Remarks</b></td></tr>" );

  line = inputFile.readline()
  while not re.search( "end OracleForOnePhaseWithShrinkingGrainSize", line):
    methodTrace        = line.split( "=")[0]
    biggestProblemSize = line.split( "=")[1].split( "," )[0]
    grainSize          = line.split( "=")[1].split( "," )[1]
    searchDelta        = line.split( "=")[1].split( "," )[3]
    accuracy           = line.split( "=")[1].split( "," )[4]
    totalSerialTime    = line.split( "=")[1].split( "," )[5]
    serialTimings      = line.split( "=")[1].split( "," )[6]
    currentTiming      = line.split( "(")[1].split( "," )[0]
    
    htmlOverview.write( "<tr>" );
    htmlOverview.write( "<td>" + methodTrace + "</td>" );
    htmlOverview.write( "<td>" + biggestProblemSize + "</td>" );
    
    analysis = Analysis.Nop
    if float(serialTimings)==0:
      analysis = Analysis.NoSerialRuntimeYet
    elif int(biggestProblemSize)==int(grainSize) and int(searchDelta)>0:
      analysis = Analysis.SeemsNotToScale
    elif int(biggestProblemSize)>int(grainSize) and int(searchDelta)>0:
      analysis = Analysis.MightScale
    elif int(biggestProblemSize)==int(grainSize) and int(searchDelta)==0:
      analysis = Analysis.DoesNotScale
    elif int(biggestProblemSize)>int(grainSize) and int(searchDelta)==0:
      analysis = Analysis.Scales


    if analysis==Analysis.NoSerialRuntimeYet:
      htmlOverview.write( "<td bgcolor=\"Yellow\">" + grainSize + "</td>" );
    elif analysis==Analysis.SeemsNotToScale:
      htmlOverview.write( "<td bgcolor=\"Fuchsia\">" + grainSize + "</td>" );
    elif analysis==Analysis.MightScale:
      htmlOverview.write( "<td bgcolor=\"LightSkyBlue\">" + grainSize + "</td>" );
    elif analysis==Analysis.DoesNotScale:
      htmlOverview.write( "<td bgcolor=\"Red\">" + grainSize + "</td>" );
    elif analysis==Analysis.Scales:
      htmlOverview.write( "<td bgcolor=\"LightGreen\">" + grainSize + "</td>" );
    else:
      htmlOverview.write( "<td bgcolor=\"White\">" + grainSize + "</td>" );

    htmlOverview.write( "<td>" + searchDelta + "</td>" );
    htmlOverview.write( "<td>" + accuracy + "</td>" );
    htmlOverview.write( "<td>" + totalSerialTime + "</td>" );
    htmlOverview.write( "<td>" + serialTimings + "</td>" );
    htmlOverview.write( "<td>" + currentTiming + "</td>" );

    if re.search( "estimated-speedup=", line):    
      htmlOverview.write( "<td>" + line.split("estimated-speedup=")[1] + "</td>" );
    elif float(currentTiming)==0:
      htmlOverview.write( "<td> </td>" );
    else:
      htmlOverview.write( "<td>not available</td>" );

    
    htmlOverview.write( "<td>" );
    if int(searchDelta)!=0:
      htmlOverview.write( "Still searching. " );
    htmlOverview.write( str(analysis) );
    htmlOverview.write( ". " );
    if analysis==Analysis.SeemsNotToScale:
      htmlOverview.write( "Code might have found scaling setup but wants to re-validate serial runtime." );
    
    costForBiggestProblemSet = 0
    if float(serialTimings)==0:
      costForBiggestProblemSet = float(currentTiming) * float(biggestProblemSize);
    if float(serialTimings)>0:
      costForBiggestProblemSet = float(serialTimings) * float(biggestProblemSize);
    htmlOverview.write( "Est. serial runtime for call on largest problem set: " + str(costForBiggestProblemSet) );
    
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

htmlOverview.write( "Number of oracles=" + str(totalNumberOfOracles) + " (incl. oracles required for repository/algorithm steering)" );
htmlOverview.write( "<br />" );
htmlOverview.write( "Oracles required for repository steering=" + str(oraclesForSteering) );

htmlOverview.write( "<h3>Table of content:</h3>" );
htmlOverview.write( "<ul>" );
for adapter in range(oraclesForSteering,totalNumberOfOracles):
  htmlOverview.write( "<li><a href=\"#adapter-" + str(adapter) + "\">Adapter " + str(adapter) + "</a></li>" );
htmlOverview.write( "</ul>" );

htmlOverview.write( "<p>Empty adapter sections imply that this adapter is not used by the code. The adapter order requals the adapter order in the specification file. Only those program phases that are actually used are also displayed.</p>" );
htmlOverview.write( "<p>All individual timings are normalised by the number of entries handled, i.e. they do specify time per grid entity. We thus may not compare them directly to the total runtime.</p>" );

for adapter in range(oraclesForSteering,totalNumberOfOracles):
  htmlOverview.write( "<h3 id=\"adapter-" + str(adapter) + "\">Adapter " + str(adapter) + "</h3>" );
  processMeasurement(adapter)
  htmlOverview.write( "<p>The serial runtime is not reliable: It tracks only timinigs if the code fragment is not parallelised and the search for a grain size still is switched on. The longer your code runs, the more `invalid' this figure becomes. However, large entries indicate that you have heavy-weight serial part in your code.</p>" );
  