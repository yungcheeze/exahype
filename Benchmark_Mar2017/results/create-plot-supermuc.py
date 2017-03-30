import sys
import re
import pylab 
import os

#import matplotlib.pyplot as pylab


Colors = [ "#ff0000", "#00ff00", "#0000ff", "#aabbaa", "#ff00ee", "#00eeff", "#ffee00" ]


MaxNodes = 8+1

CoresPerNode = 28

def plot(ranksPerNode,variant):
 filename = sys.argv[1] + "-" + str(ranksPerNode) + "-ranks-per-node-" + variant
 print filename

 Nodes = [ 1, 2, 4, 6, 8 ] 
 
 pylab.clf()
 serialRuntimeMax = 0
 serialRuntimeMin = 65536
 nodesMax         = 0

 Ps         = [3,9]
 Extensions = ["no-output-notbb", "output-notbb", "no-output-tbb", "output-tbb", "no-output-mtmpi", "output-mtmpi"]

 foundData = False

 for p in Ps:
  for extension in Extensions:
   xData = []
   yData = []  
   for nodes in range(0,MaxNodes):
    inputFilename = sys.argv[1] + "-" + str(ranksPerNode) + "-ranks-per-node-" + str(nodes) + "-nodes-" + variant + "-p" + str(p) + "-" + extension + ".out"
    if os.path.isfile( inputFilename ):
     print "read " + inputFilename
     foundData = True
     try:
      inputFile = open( inputFilename,  "r" )
      for line in inputFile:
       if re.search( "ADERDGTimeStep ", line ):
        time       = float(line.strip().split( "|" )[6])
        if nodes==1:
         xData.append( 1 )
        else:
         xData.append( nodes*CoresPerNode )
        yData.append( time )
     except Exception as e:
      print "failure while parsing " + inputFilename + ": " + str(e)  

    if len(yData)>0 and serialRuntimeMax<yData[0] and xData[0]==1:
      serialRuntimeMax = yData[0]
    if len(yData)>0 and serialRuntimeMin>yData[0] and xData[0]==1:
      serialRuntimeMin = yData[0]
    if len(xData)>0 and nodesMax<xData[-1]:
       nodesMax = xData[-1]  

   Colors  = [ "#ff0000", "#00ff00", "#0000ff", "#ffff00", "#ff00ff", "#00ffff" ]
   Symbols = [ "o",       "s",       "<",       ">",       "^",       "v" ] 
   myMarkerSize = 8
   myAlpha      = 0.5
   myColor      = Colors[ Extensions.index(extension) ]
   
   mySymbol     = ""
   if Ps.index(p)==0:
     mySymbol   = "-"
   else:
     mySymbol   = "--"
   mySymbol     = mySymbol + Symbols[ Extensions.index(extension) ]  

   if len(xData)>0 and re.search( "no-output", extension):
     pylab.plot( xData, yData, mySymbol, label=extension + ", p=" + str(p), color=myColor, markerfacecolor='none', markersize=myMarkerSize, alpha=myAlpha)
   elif len(xData)>0:
     pylab.plot( xData, yData, mySymbol, label=extension + ", p=" + str(p), color=myColor,                         markersize=myMarkerSize, alpha=myAlpha)

 if nodesMax>0 and serialRuntimeMin<65536:
  xData = [1,nodesMax]
  yData = [serialRuntimeMin,serialRuntimeMin/(nodesMax*CoresPerNode)]
  pylab.plot( xData, yData, "--", label="serial", color="#787878")

 if nodesMax>0 and serialRuntimeMax>0:
  xData = [1,nodesMax]
  yData = [serialRuntimeMax,serialRuntimeMax/(nodesMax*CoresPerNode)]
  pylab.plot( xData, yData, "--", color="#787878")
  
 pylab.title( sys.argv[1] + ", experiment " + variant + ", " + str(ranksPerNode) + " ranks/node" )
 pylab.xlabel('cores')
 pylab.ylabel('runtime/grid sweep [t]=s')
 


 if foundData:
  pylab.legend(fontsize=10, ncol=2, framealpha=0.5, loc='best')
  pylab.savefig( filename + ".png")
  pylab.savefig( filename + ".pdf")
  try:
    pylab.legend().set_visible(False)
  except:
    pass
  pylab.savefig( filename + "-no-legend.png")
  pylab.savefig( filename + "-no-legend.pdf")
  pylab.legend(fontsize=10, ncol=2, framealpha=0.5, loc='best')
  pylab.yscale( 'log', basey=10 )  
  pylab.savefig( filename + "-log.png")
  pylab.savefig( filename + "-log.pdf")
  try:
    pylab.legend().set_visible(False)
  except:
    pass
  pylab.savefig( filename + "-log-no-legend.png")
  pylab.savefig( filename + "-log-no-legend.pdf")


if len(sys.argv)!=2:
  print "please specify file prefix (such as Z4)"
  exit(-1)


for ranksPerNode in [ 1, 2, 4, 7, 14, 28 ]:
  for variant in [ "0", "1", "2", "3", "4" ]:
    plot(ranksPerNode,variant)
