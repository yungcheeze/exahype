import sys
import re
import pylab 
import os




def processMeasurement(adapter,phase):
  substring    = line.split( "adapter=" + str(adapter) + ", phase=" + str(phase) )[1]
  substring    = substring.split( "adapter" )[0]

  name         = substring.split("name=")[1].split(",")[0] 
  problemSize  = substring.split("problem-size=")[1].split(")")[0] 
  grainSizes   = int(substring.split("[")[1].split( "grain size(s) studied" )[0])

  htmlOverview.write( "name=" + name + "<br />" );
  htmlOverview.write( "max problem size=" + problemSize + "<br />" );
  htmlOverview.write( "grain sizes studied=" + str(grainSizes) + "<br />" );

  
  xValues = []
  yValues = []
  zValues = []
  
  for measurement in substring.split( "(grain-size=" )[1:-1]:
    grainSize  = float(measurement.split( "," )[0])
    time       = float(measurement.split( "(" )[1].split( "," )[0])
    deviation  = float(measurement.split( "deviation=" )[1].split( ")" )[0])
    if xValues==[] or xValues[-1]<grainSize:
      #htmlOverview.write( "grain size=" + str(grainSize) + "<br />" );
      #htmlOverview.write( "time=" + str(time) + "<br />" );
      #htmlOverview.write( "deviation=" + str(deviation) + "<br />" );
      xValues.append( grainSize )
      yValues.append( time )
      zValues.append( deviation  )

  pylab.clf()
  outputFileName = "adapter-" + str(adapter) + "-phase-" + str(phase)
  print "write file " + outputFileName
    
  symbolAndColourCounter = 0
  
  try:     
    #, alpha=AlphaValue, label=solver
    pylab.plot(xValues, yValues, "-s",  markersize=10, color="#0000ff" )
    pylab.plot(xValues, zValues, "o",   markersize=12, color="#ff0000" )
    pylab.grid(True)
    #pylab.title( "$h_{max}$=" + "%-4.2e" % hMax + ", $h_{min}$=" + "%-4.2e" % hMin + ", $\omega=$" + str(omega) + ", $\\theta=$" + str(theta) )
    pylab.ylabel('[t]=s')
    pylab.yscale( 'log' )
    #pylab.ylim([1e-14,2])
    pylab.xlabel('grain size')
    pylab.legend(fontsize=9, framealpha=0.5)
    pylab.legend(loc='upper left',framealpha=0.5)
  
    pylab.savefig( outputFileName + ".png")
    pylab.savefig( outputFileName + ".pdf")
    htmlOverview.write( "<img src=\"" + outputFileName + ".png\" /> <br />" );
  except:
    print "Unexpected error:", sys.exc_info()[0]







  #  symbolAndColourCounter = symbolAndColourCounter + 1
    




#
#
#   main
# ========
#
#    
if (len(sys.argv)!=4):
  print "usage: python ../postprocess-sampling-output.py outputfile no-of-adapter no-of-phases"
  print "  no-of-adapter  number of adapters you have in your project (see your specification file)"
  print "  no-of-phases   number of code phases that are tuned via an oracle. 19 by default (cmp OracleForOnePhase)"
  quit()

htmlOverview = open( sys.argv[1] + ".html",  "w" )
htmlOverview.write( "<h1>" + sys.argv[1] + "</h1>" );

inputFile = open(sys.argv[1], "r" )
line      = inputFile.readline()
numberOfAdapters = int(sys.argv[2])
numberOfPhases   = int(sys.argv[3])

htmlOverview.write( "Number of adapters=" + str(numberOfAdapters) );
htmlOverview.write( "<br />" );
htmlOverview.write( "Number of phases=" + str(numberOfPhases) );

htmlOverview.write( "<p>The blue lines present timings for particular grain sizes. The red dots denote the standard deviation belonging to the measurements.</p>" );

for adapter in range(0,numberOfAdapters):
  htmlOverview.write( "<h2>Adapter " + str(adapter) + "</h2>" );
  for phase in range(0,numberOfPhases):
    htmlOverview.write( "<h3>Phase " + str(phase) + "</h3>" );
    try:
      processMeasurement(adapter,phase)
    except:
      pass
