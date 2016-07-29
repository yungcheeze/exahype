import argparse
from argparse import RawTextHelpFormatter

import pylab
import runtimeParser



########################################################################
# START OF THE PROGRAM
########################################################################
# For the meaning of the command line arguments, see the documentation
# of the function 'plot_multithreading_adapter_scaling' above.
help = '''
Creates a speedup plot based on Peano output files with specific file naming pattern.
If multiple adapters are specified, then the cumulative user times are used to compute the speedup.
\n\n
Sample usage:\n
python usertimeplot.py -path \'examples/151217_phi1_node\' -prefix \'151217\' -legend \'2x Xeon  5-2650 @ 2.00GHz\' -mode TBB -cc icpc -ylim 16 -per_iteration -adapter \'Predictor+Corrector\' -t 1 2 4 6 8 10 12 16 32 -hyperthreading'''

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('-table',nargs='+',required=True,help="Tables that are to be read.")
parser.add_argument('-output',required=True,help="Output file (extensions pdf and png are added).")
parser.add_argument('-adapter',default=['Total'],nargs='+',help="Name of the adapters separated by a \'+\'. (Use \'Total\' for the  cumulative time of all specified adapters. Does not make sense with \'per_iteration\' switched on.)")
args   = parser.parse_args();


# @todo Dimension in **2


pylab.clf()
pylab.xlabel('cores')


for table in args.table:
  print "read " + table
  xdata    = runtimeParser.readColumnFromTable(table,0)
  maxLevel = runtimeParser.readColumnFromTable(table,1)   

  for adap in args.adapter:
    totalTime    = runtimeParser.readColumnFromTable(table, runtimeParser.getAdapterRuntimeColumnFromTable(table,adap) )
    count        = runtimeParser.readColumnFromTable(table, runtimeParser.getAdapterCountColumnFromTable(table,adap) )
    ydata = []
    for i in range(0,len(totalTime)):
      ydata.append( totalTime[i]/(3**(maxLevel[i]*2))/count[i])
    pylab.plot(xdata,ydata,markersize=4,label=table,marker='',markevery=1,lw=1.2,linestyle='dashed',color='grey') 

pylab.ylabel( "time [t]=s (normalised by finest regular grid)" )

pylab.savefig( args.output + ".png" )
pylab.savefig( args.output + ".pdf" )
#pylab.savefig( args.output + ".eps" )

