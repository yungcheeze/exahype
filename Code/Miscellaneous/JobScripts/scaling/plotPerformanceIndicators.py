import matplotlib
matplotlib.use('Agg') 

import pylab
import sys

def loadData():
  filename = sys.argv[1]
  file = open(filename, "r" )
    
  line = file.readline()
  assert(line == "\"JobId\",\"Application\",\"Dim\",\"nProcs\",\"p\",\"hMax\",\"Compiler\",\"Mode\",\"optimisticADERDGSteps\",\"averageTimestepCostOfSingleADERDGStep\",\"numberOfSteps\",\"simulationRuntime\",\"costOfADERDGSteps\"\n")
  
  line = line.replace("\"", "")
  
  
  columns = line.split(",")
  print columns
  groupByColumn = [1, 2, 4, 5, 6, 7]
  plots = []
  speedups = []
  
  data = [];
  for line in file:
    newRow = line.split(",")
    
    if len(newRow) != 13:
      continue
      
      
    # print newRow
    if newRow[9] == "-1":
      continue
      
    
    data.append(newRow)
  
  
    plotName = ""
    for column in groupByColumn:
      plotName += columns[column] + "=" + str(newRow[column]) + ";"

    if plotName not in plots:
      plots.append(plotName)
      speedups.append([])
      print plots.index(plotName)
      
    speedups[plots.index(plotName)].append([newRow[3], newRow[8], newRow[9], newRow[10], newRow[11], newRow[12]])
        
        
      
      
  for plotName in plots:
    writeOutPlot(plotName, speedups[plots.index(plotName)])

def writeOutPlot(plotName, measurements):

  pylab.clf()

  experimentSetCounter = 5
  Colors           = ['#ff0000','#00ff00','#0000ff','#ffff00','#ff00ff','#00ffff','#00ffff','#00ffff','#00ffff','#00ffff','#00ffff','#00ffff']
  Markerfacecolors = ['None','None','None','k','k','k','k','k','k','k','k','k','k']
  Markers          = ['o','^','s','o','^','s','s','s','s','s','s','s','s']
  adap = "test"

  print plotName
  print measurements
  symbolCounter = 0
  scalingFactor = -1
  
  for measurement in measurements:
    print "measurement[0] ", measurement[0]
    if measurement[0] == "2":
      scalingFactor = float(measurement[2])
    
  
  print scalingFactor
  
  for measurement in measurements:
    symbolCounter = 1
    # print measurement
    print measurement[0], measurement[2], abs(float(measurement[2]))
    pylab.plot(measurement[0], float(measurement[2]), markersize=experimentSetCounter+4,label=adap,color=Colors[symbolCounter],marker=Markers[symbolCounter],markerfacecolor=Markerfacecolors[symbolCounter],markevery=1,lw=1.2) 

  pylab.loglog( basex=2, basey=2 )

  
  pylab.ylabel( "time [t]=s" )
  pylab.xlabel( "Kekse" )
  pylab.grid(True)
  pylab.xlabel('cores')
  pylab.savefig(plotName + ".runtime.png")

  symbolCounter = 0
  pylab.clf()
  if scalingFactor != -1:
    for measurement in measurements:
      symbolCounter = 1
      # print measurement
      print measurement[0], float(measurement[2])/scalingFactor
      if (measurement[0] != "1"):
        pylab.plot(measurement[0], float(measurement[2])/(scalingFactor*(float(measurement[0])-1)), markersize=experimentSetCounter+4,label=adap,color=Colors[symbolCounter],marker=Markers[symbolCounter],markerfacecolor=Markerfacecolors[symbolCounter],markevery=1,lw=1.2) 
      else:
        pylab.plot(measurement[0], float(measurement[2])/(scalingFactor), markersize=experimentSetCounter+4,label=adap,color=Colors[symbolCounter],marker=Markers[symbolCounter],markerfacecolor=Markerfacecolors[symbolCounter],markevery=1,lw=1.2) 

    pylab.loglog( basex=2, basey=2 )

    
    pylab.ylabel( "time [t]=s" )
    pylab.xlabel( "Kekse" )
    pylab.grid(True)
    pylab.xlabel('cores')
    pylab.savefig(plotName + ".efficiency.png")

  print
  print
  print
  print
  
  
loadData()
