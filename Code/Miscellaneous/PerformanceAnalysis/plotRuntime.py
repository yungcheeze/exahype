import argparse
from argparse import RawTextHelpFormatter

import pylab
import runtimeParser


Colors           = ['#ff0000','#00ff00','#0000ff','#ffff00','#ff00ff','#00ffff']
Markerfacecolors = ['None','None','None','k','k','k']
Markers          = ['o','^','s','o','^','s']


xDataMax = 0


def addData(table,normalisation,plotLabels,experimentSetCounter,label):
  global xDataMax
  
  xdata    = runtimeParser.readColumnFromTable(table,0)

  symbolCounter = 0
  yDataMin      = 65636
  for adap in args.adapter:
    totalTime    = runtimeParser.readColumnFromTable(table, runtimeParser.getAdapterRuntimeColumnFromTable(table,adap) )
    count        = runtimeParser.readColumnFromTable(table, runtimeParser.getAdapterCountColumnFromTable(table,adap) )
    ydata = []
    for i in range(0,len(totalTime)):
      if count[i]==0:
        print "ERROR " + str(i) + "th entry in " + table + "'s " + str(runtimeParser.getAdapterCountColumnFromTable(table,adap)) + "th column (adapter " + adap + ") equals 0"
        ydata = []
        break 
      ydata.append( totalTime[i]*normalisation/count[i])
    if len(ydata)==0:
      print "WARNING: file " + table + " seems to be empty for adapter " + adap
    else:
      if ydata[0]<yDataMin:
        yDataMin = ydata[0]
      if xdata[-1]>xDataMax:
        xDataMax = xdata[-1]
      if (plotLabels):
        pylab.plot(xdata,ydata,markersize=experimentSetCounter+4,label=adap,color=Colors[symbolCounter],marker=Markers[symbolCounter],markerfacecolor=Markerfacecolors[symbolCounter],markevery=1,lw=1.2) 
      else:
        pylab.plot(xdata,ydata,markersize=experimentSetCounter+4,color=Colors[symbolCounter],marker=Markers[symbolCounter],markerfacecolor=Markerfacecolors[symbolCounter],markevery=1,lw=1.2) 
    symbolCounter = symbolCounter + 1

  if len(xdata)>0 and len(ydata)>0:
    pylab.text(xdata[-1],ydata[-1],label)

  ydata = [yDataMin/x*xdata[0] for x in xdata]
  pylab.plot(xdata,ydata,markersize=4,markevery=1,lw=1.2,linestyle='dashed',color='grey') 



  
def initGlobalPlotterSettings():
  pylab.ylabel( "time [t]=s" )
  pylab.xlabel( args.xaxislabel )
  pylab.grid(True)
  pylab.legend(loc='best',fontsize='%d' % int(args.fontsize))  
  #fig, ax = pylab.subplots()
  #ax.yaxis.set_major_formatter(pylab.FormatStrFormatter('%.2e'))

  pylab.ticklabel_format(axis='y',style='sci', scilimits=(-2,2))
  #ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))


def switchToLogScales():
  global xDataMax
  
  pylab.loglog( basex=2, basey=2 )
  XTicks  = [1]
  XLabels = [ "serial" ]
  for i in range(1,int(xDataMax)+2):
    if i != 0 and ((i & (i - 1)) == 0):
      XTicks.append( i )
      XLabels.append( str(i) ) 
  pylab.xticks(XTicks,XLabels)


########################################################################
# START OF THE PROGRAM
########################################################################
# For the meaning of the command line arguments, see the documentation
# of the function 'plot_multithreading_adapter_scaling' above.
help = '''
Creates a speedup plot based on Peano output files with specific file naming pattern.
If multiple adapters are specified, then the cumulative user times are used to compute the speedup.
\n\n
Sample usages:\n
python ../plotRuntime.py -experimentdescription '' -table x4-1.results.table -output experiment1 -xaxislabel "MPI Ranks" -dimension 2 \n
python ../plotRuntime.py -experimentdescription '' -table x4-1.results.table -adapter Total ADERDGTimeStep PredictorRerun -output experiment1 -xaxislabel "MPI Ranks" -dimension 2 \n
python ../plotRuntime.py -experimentdescription 'depth 6' 'depth 7' 'depth 8' -table x4-1.results.table x4-2.results.table x4-3.results.table -adapter ADERDGTimeStep PredictorRerun -output experimentx4 -xaxislabel "MPI Ranks" -dimension 2 \n
python ../plotRuntime.py -experimentdescription 'depth 4' 'depth 5' 'depth 6' 'depth 7' 'depth 8' 'depth 9' -table x16-0.results.table x16-1.results.table x16-2.results.table x16-3.results.table x16-4.results.table x16-5.results.table -adapter ADERDGTimeStep PredictorRerun -output experimentx16 -xaxislabel "MPI Ranks" -dimension 2 \n
'''

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('-table',nargs='+',required=True,help="Tables that are to be read.")
parser.add_argument('-output',required=True,help="Output file (extensions pdf and png are added).")
parser.add_argument('-adapter',default=['Total'],nargs='+',help="Name of the adapters separated by a \'+\'. (Use \'Total\' for the  cumulative time of all specified adapters. Does not make sense with \'per_iteration\' switched on.)")
parser.add_argument('-dimension',required=True,help="Dimension of problem. Either 2 or 3.")
parser.add_argument('-xaxislabel',required=True,help="Label of x axis.")
parser.add_argument('-fontsize',default=10,required=False,help="Font size of the legend and tick labels. Axis labels are computed by ceiling the font size times a factor 1.2.")
parser.add_argument('-experimentdescription',nargs='+',required=True,help="Per table entry, one experiment desription is required")
args   = parser.parse_args();

dim = int(args.dimension)



#
# Raw runtime
#
outputFile = args.output + "-raw-runtime"
pylab.clf()
pylab.xlabel('cores')

experimentSetCounter =  0
for (table,label) in zip(args.table,args.experimentdescription):
  print "read " + table
  maxLevel = runtimeParser.readColumnFromTable(table,1)   
  addData(table,1.0,table==args.table[0],experimentSetCounter,label)
  experimentSetCounter = experimentSetCounter + 1

initGlobalPlotterSettings()

pylab.savefig( outputFile + ".png" )
pylab.savefig( outputFile + ".pdf" )
switchToLogScales()
pylab.savefig( outputFile + "-log.png" )
pylab.savefig( outputFile + "-log.pdf" )
print "written " + outputFile




#
# Runtime scaled by regular grid of max depth
#
outputFile = args.output + "-scaled-by-regular-grid"
pylab.clf()
pylab.xlabel('cores')

experimentSetCounter =  0
for (table,label) in zip(args.table,args.experimentdescription):
  print "read " + table
  maxLevel = runtimeParser.readColumnFromTable(table,1)   
  if len(maxLevel)>0:
    normalisation = 1.0/(3**(maxLevel[-1]*dim))
    addData(table,normalisation,table==args.table[0],experimentSetCounter,label)
  else:
    print "WARNING: Could not determine max level for table " + table
  experimentSetCounter = experimentSetCounter + 1

initGlobalPlotterSettings()

pylab.savefig( outputFile + ".png" )
pylab.savefig( outputFile + ".pdf" )
switchToLogScales()
pylab.savefig( outputFile + "-log.png" )
pylab.savefig( outputFile + "-log.pdf" )
print "written " + outputFile
    