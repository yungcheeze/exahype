import sys
import re
import pylab 
import os



Colours   = [ "#ff0000", "#00ff00", "#0000ff", "#ffff00", "#ff00ff", "#00ffff", "#232323", "#aaaaaa" ]
Symbols   = [ "-s",      "-o",      "-^",      "-v",      "-<",      "->",      "-D",      "-1" ]



def processMeasurement(adapter,phase):
  substring    = line.split( "adapter=" + str(adapter) + ", phase=" + str(phase) )[1]

  name         = substring.split("name=")[1].split(",")[0] 
  problemSize  = substring.split("problem-size=")[1].split(")")[0] 
  grainSizes   = int(substring.split("[")[1].split( "grain size(s) studied" )[0])

  htmlOverview.write( "name=" + name + "<br />" );
  htmlOverview.write( "max problem size=" + problemSize + "<br />" );
  htmlOverview.write( "grain sizes studied=" + str(grainSizes) + "<br />" );
  
  for measurement in substring.split( "(grain-size=" )[1:-1]:
    #htmlOverview.write( "<pre>" );
    #htmlOverview.write( measurement  );
    #htmlOverview.write( "</pre>" );
    grainSize  = float(measurement.split( "," )[0])
    htmlOverview.write( "grain size=" + str(grainSize) + "<br />" );
    time       = float(measurement.split( "(" )[1].split( "," )[0])
    htmlOverview.write( "time=" + str(time) + "<br />" );
    deviation  = float(measurement.split( "deviation=" )[1].split( ")" )[0])
    htmlOverview.write( "deviation=" + str(deviation) + "<br />" );
   


#
#
#   Compare the different solvers
# =================================
#
#    
def plotSolverDependencies(hMax,hMin,omega,theta):
  pylab.clf()
  
  outputFileName = "solvers-" + str(hMax) + "-" + str(hMin) + "-" + str(omega) + "-" + str(theta) + "-res2"
  
  print "create " + outputFileName 

  symbolAndColourCounter = 0
  
  for solver in Solvers:
    inputFile = getFilename(hMax,hMin,solver,omega,theta) + ".table"
    pylab.plot(readColumnFromTable(inputFile,0), readColumnFromTable(inputFile,5), Symbols[symbolAndColourCounter], markersize=10, label=solver, color=Colours[symbolAndColourCounter], markevery=MarkEveryOffset-symbolAndColourCounter, alpha=AlphaValue)
    symbolAndColourCounter = symbolAndColourCounter + 1
    
  try:
    pylab.grid(True)
    pylab.title( "$h_{max}$=" + "%-4.2e" % hMax + ", $h_{min}$=" + "%-4.2e" % hMin + ", $\omega=$" + str(omega) + ", $\\theta=$" + str(theta) )
    pylab.ylabel('$|r(n)|_{2}/|r(0)|_{2}$')
    pylab.yscale( 'log' )
    pylab.ylim([1e-14,2])
    pylab.xlabel('iteration $n$')
    pylab.legend(fontsize=9, framealpha=0.5)
    pylab.legend(loc='lower left',framealpha=0.5)
  
    pylab.savefig( convertIntoImageFileName(outputFileName) + ".png")
    pylab.savefig( convertIntoImageFileName(outputFileName) + ".pdf")
    htmlOverview.write( "<img src=\"" + convertIntoImageFileName(outputFileName) + ".png\" />" );
    pylab.legend().set_visible(False)
    pylab.savefig( convertIntoImageFileName(outputFileName) + "-no-legend.png")
    pylab.savefig( convertIntoImageFileName(outputFileName) + "-no-legend.pdf")
  except:
    print "Unexpected error:", sys.exc_info()[0]
        

  pylab.clf()
  
  outputFileName = "solvers-" + str(hMax) + "-" + str(hMin) + "-" + str(omega) + "-" + str(theta) + "-resmax"
  
  print "create " + outputFileName 

  symbolAndColourCounter = 0
  
  for solver in Solvers:
    inputFile = getFilename(hMax,hMin,solver,omega,theta) + ".table"
    pylab.plot(readColumnFromTable(inputFile,0), readColumnFromTable(inputFile,6), Symbols[symbolAndColourCounter], markersize=10, label=solver, color=Colours[symbolAndColourCounter], markevery=MarkEveryOffset-symbolAndColourCounter, alpha=AlphaValue)
    symbolAndColourCounter = symbolAndColourCounter + 1

  try:     
    pylab.grid(True)
    pylab.title( "$h_{max}$=" + "%-4.2e" % hMax + ", $h_{min}$=" + "%-4.2e" % hMin + ", $\omega=$" + str(omega) + ", $\\theta=$" + str(theta) )
    pylab.ylabel('$|r(n)|_{max}/|r(0)|_{max}$')
    pylab.yscale( 'log' )
    pylab.ylim([1e-14,2])
    pylab.xlabel('iteration $n$')
    pylab.legend(fontsize=9, framealpha=0.5)
    pylab.legend(loc='lower left',framealpha=0.5)
  
    pylab.savefig( convertIntoImageFileName(outputFileName) + ".png")
    pylab.savefig( convertIntoImageFileName(outputFileName) + ".pdf")
    htmlOverview.write( "<img src=\"" + convertIntoImageFileName(outputFileName) + ".png\" /> <br />" );
    htmlOverview.write( "<a href=\"#toc\">Up to table of content</a><br />" );
    pylab.legend().set_visible(False)
    pylab.savefig( convertIntoImageFileName(outputFileName) + "-no-legend.png")
    pylab.savefig( convertIntoImageFileName(outputFileName) + "-no-legend.pdf")
  except:
    print "Unexpected error:", sys.exc_info()[0]



#
#
#   main
# ========
#
#    
if (len(sys.argv)!=2):
  print "usage: python ../postprocess-sampling-output.py outputfile"
  quit()

htmlOverview = open( sys.argv[1] + ".html",  "w" )
htmlOverview.write( "<h1>" + sys.argv[1] + "</h1>" );

inputFile = open(sys.argv[1], "r" )
line      = inputFile.readline()
numberOfAdapters = int(line.split( "adapters:" )[1].split(",")[0])
numberOfMethods  = int(line.split( "phases:" )[1].split(",")[0])

htmlOverview.write( "Number of adapters=" + str(numberOfAdapters) );
htmlOverview.write( "<br />" );
htmlOverview.write( "Number of phases=" + str(numberOfMethods) );

for adapter in range(0,numberOfAdapters):
  htmlOverview.write( "<h2>Adapter " + str(adapter) + "</h2>" );
  for phase in range(0,numberOfMethods):
    htmlOverview.write( "<h3>Method " + str(phase) + "</h3>" );
    processMeasurement(adapter,phase)
#  processAdapter()
#for hMin in Hs:
#  for omega in Omegas:
#    for theta in Thetas:
#      hMax = Hs[0]
#      currentHMin = hMin
#      if (sys.argv[2]=="only-adaptive-mg"):
#        currentHMin = hMax
#      while hMax >= currentHMin:
#        plotSolverDependencies(hMax,hMin,omega,theta)
#        hMax = hMax / 3