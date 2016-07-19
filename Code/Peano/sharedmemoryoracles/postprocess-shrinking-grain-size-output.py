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
  
  htmlOverview.write( "<h3 id=\"adapter-" + str(adapter) + "\">Adapter " + str(adapter) + "</h3>" );
  htmlOverview.write( "<table border=\"1\">" );
  htmlOverview.write( "<tr><td><b>Method</b></td><td><b>Grain Size</b></td><td><b>Maximum problem size</b></td><td><b>Standard deviation</b></td><td><b>Remarks</b></td></tr>" );

  inputFile = open(sys.argv[1], "r" )
  for line in inputFile:
    try:
      isRightAdapter = re.search(searchPattern,line)
      if isRightAdapter and re.search("still determining serial runtime",line):
        htmlOverview.write( "<tr>" );
        htmlOverview.write( "<td>" + line.split( "method=")[1].split(":")[0] + "</td>" );
        htmlOverview.write( "<td bgcolor=\"red\">still searching for serial runtime, no meaningful data available</td>" );
        htmlOverview.write( "<td />" );
        htmlOverview.write( "<td>" + line.split( "std-deviation=")[1].split(")")[0] + "</td>" );
        htmlOverview.write( "<td>Accurate value is value where new measurements do not change std deviation anymore by more  " + line.split( "eps=")[1].split(",")[0] + "</td>" );
        htmlOverview.write( "</tr>" );
      if isRightAdapter and re.search("does not scale, oracle is not searching anymore",line):
        htmlOverview.write( "<tr>" );
        htmlOverview.write( "<td>" + line.split( "method=")[1].split(":")[0] + "</td>" );
        htmlOverview.write( "<td bgcolor=\"white\">does not scale, oracle is not searching anymore</td>" );
        htmlOverview.write( "<td />" );
        htmlOverview.write( "<td>" + line.split( "std-deviation=")[1].split(")")[0] + "</td>" );
        htmlOverview.write( "<td>Accurate value is value where new measurements do not change std deviation anymore by more  " + line.split( "eps=")[1].split(",")[0] + "</td>" );
        htmlOverview.write( "</tr>" );
      if isRightAdapter and re.search("still searching for optimal grain size",line):
        htmlOverview.write( "<tr>" );
        htmlOverview.write( "<td>" + line.split( "method=")[1].split(":")[0] + "</td>" );
        htmlOverview.write( "<td bgcolor=\"yellow\">searching (current grain size " );
        htmlOverview.write( line.split( "currentGrainSize=" )[1].split(",")[0] + "</td>" );
        htmlOverview.write( "<td>" + line.split( "biggestProblemSize=" )[1].split(",")[0] + "</td>" );
        htmlOverview.write( "<td>" + line.split( "std-deviation=")[1].split(")")[0] + "</td>" );
        htmlOverview.write( "</tr>" );
      if isRightAdapter and re.search("scales, oracle is not searching anymore",line):
        htmlOverview.write( "<tr>" );
        htmlOverview.write( "<td>" + line.split( "method=")[1].split(":")[0] + "</td>" );
        htmlOverview.write( "<td bgcolor=\"green\">does scale</td>" );
        htmlOverview.write( "<td>" + line.split( "currentGrainSize=" )[1].split(",")[0] + "</td>" );
        htmlOverview.write( "<td>" + line.split( "biggestProblemSize=" )[1].split(",")[0] + "</td>" );
        htmlOverview.write( "<td>" + line.split( "std-deviation=")[1].split(")")[0] + "</td>" );
        htmlOverview.write( "</tr>" );
        
    except:
      print "ERROR processing adapter " + str(adapter)

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
numberOfAdapters      = int( line.split("=")[1] )
line                  = inputFile.readline()
adaptersForSteering   = int( line.split("=")[1] )
line                  = inputFile.readline()
noOfMethodsCalling    = int( line.split("=")[1] )
inputFile.close()

htmlOverview.write( "Total number of oracles=" + str(totalNumberOfOracles) );
htmlOverview.write( "<br />" );
htmlOverview.write( "Number of adapters=" + str(numberOfAdapters) + " (incl. adapters required for algorithm steering)" );
htmlOverview.write( "<br />" );
htmlOverview.write( "Adapters required for repository steering=" + str(adaptersForSteering) );
htmlOverview.write( "<br />" );
htmlOverview.write( "Number of methods calling=" + str(noOfMethodsCalling) + " (max)");

htmlOverview.write( "<h3>Table of content:</h3>" );
htmlOverview.write( "<ul>" );
for adapter in range(0,numberOfAdapters-adaptersForSteering):
  htmlOverview.write( "<li><a href=\"#adapter-" + str(adapter) + "\">Adapter " + str(adapter) + "</a></li>" );
htmlOverview.write( "</ul>" );

for adapter in range(0,numberOfAdapters-adaptersForSteering):
  processMeasurement(adapter)
