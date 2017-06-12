#!/usr/bin/python3
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
from argparse import RawTextHelpFormatter

import matplotlib.pylab
import matplotlib.pyplot as plt

import operator

import csv

import runtimeParser as rp
import hpclib as hpc
from plotting import scalingplot as sp

def query_table(filename,adapter,cores): 
    '''
    Read a table and filter out certain
    columns 
    Args:
    selector
       a lambda expression
    '''
    datafile    = open(filename, 'r')
    reader      = csv.reader(datafile,delimiter='&')
    sorted_data = list(filter(lambda x : x[3]==adapter and x[2] in cores, reader))
    datafile.close() 
    return sorted_data


'''
.. module:: usertimeplot
  :platform: Unix, Windows, Mac
  :synopsis: Creates a speedup plot based on Peano output files with specific file naming pattern.
   
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>

:synopsis: Creates a speedup plot based on user and CPU times stored in CSV tables.
'''

def plot_multithreading_adapter_scaling(tables,prefixes,legends,adapters,cores,n_runs,xticks,yticks,ylim,per_iteration=False,hyperthreading=False,annotate=False,create_plot=False,_fontsize=10):
    '''
    Creates a scaling plot for the cumulative user time spent within the 
    specified adapters.
   
    Args:
      tables (str[]):
         A bunch of csv files
      prefixes (str[]):
         Prefix of the files - usually the date of the test and an identifier for the machine and the MPI process that has written the output files. Must be supplied once per 'root_dir' entry.
      legends (str[]):
         Legend entry for the each data set - usually an identifier for the machine and the MPI process that has written the output files. Must be supplied once per 'root_dir' entry.
      adapters (str[]):
         Name of the adapters. (Use 'Total' for the  cumulative time forall adapters. Does not make sense with per_iteration switched on.)
      cores (int[]):
         The CPU cores (threads) you want to plot.
      xticks (str[]):
        The x-ticks.
      yticks {str[]):
        The y-ticks.
      ylim (float):
         Upper limit for the y-Axis.
      hyperthreading (bool):
         The last thread count corresponds to a hyperthreading run.
      annotate (bool):
         Annotate the plots with the speedup values.
      per_iteration (bool):
         Use the adapter times per iteration.
    ''' 
    colors           = ['k','k','k','k','k','k']
    markerfacecolors = ['None','None','None','k','k','k']
    markers          = ['o','^','s','o','^','s']
    
    n_tables = len(tables)
    n_cores  = len(cores)
    
    # Plotting
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    
    p       = list(map(int,cores))
    p2      = p
    p_ticks = list(map(str,p))

    p_ideal = [float(x) / float(p[0]) for x in p]
    
    # Ideal speedup
    # In the hyper-threading case, we only want to plot a line
    # for the 'real' threads.
    if not hyperthreading: 
        plt.plot(p,p_ideal,label=r'ideal',markersize=4,marker='',markevery=1,lw=1.2,linestyle='dashed',color='grey')
    else:
        p2        = p[0:-1] # [inclusive:exclusive]
        p2.append(p2[-1]+2)
        p_ticks[-1] = '%d+HT' % p[-2]
        plt.plot(p[0:-1],p_ideal[0:-1],label=r'ideal',markersize=4,marker='',markevery=1,lw=1.2,linestyle='dashed',color='grey')
    
    t_ref = 0
    
    # Loop over tables
    for i in range(0,n_tables):
        cumulative_user_times = [0.0]*n_cores
        cumulative_cpu_times  = [0.0]*n_cores
        iterations = 1.0
        for adapter in adapters[i]:
            data = query_table(tables[i],adapter,cores)
            iterations            = float(data[4][0]) # must be the same for all adapters
            cumulative_user_times = list(map(lambda x,y : x+y, cumulative_user_times, list(map(float,[row[5] for row in data]))))
            cumulative_cpu_times  = list(map(lambda x,y : x+y, cumulative_cpu_times,  list(map(float,[row[6] for row in data]))))

        if periteration:
            cumulative_user_times[:] = [x / iterations for x in cumulative_user_times]
            cumulative_cpu_times[:]  = [x / iterations for x in cumulative_cpu_times]
        print(cumulative_user_times)
        print(cumulative_cpu_times)
 
        if i==0:
            t_ref = cumulative_user_times[0]
        
        speedup_measured = hpc.compute_speedup(cumulative_user_times,t_ref)
    
        # Measured speedup
        sp.plot_scaling(ax,p,speedup_measured,legends[i],colors[i],markers[i],markerfacecolors[i],hyperthreading,annotate)

    plt.ylabel(r'speedup',         fontsize=float(1.2*float(_fontsize)))    
    plt.xlabel(r'number of cores', fontsize=float(1.2*float(_fontsize)))
    plt.grid(True)

    p_ticks_filtered = p_ticks
    for ip in range(0,len(p_ticks)):
        if 'HT' not in p_ticks_filtered[ip]:
            if p_ticks_filtered[ip] not in xticks:
                p_ticks_filtered[ip] = ''
    plt.xticks(p2,p_ticks_filtered)

    p_ideal_filtered = list(map(int,p_ideal))
    p_ideal_filtered = list(map(str,p_ideal_filtered))
    for ip in range(0,len(p_ideal)):
       if p_ideal_filtered[ip] not in yticks:
         p_ideal_filtered[ip] = ''
    plt.yticks(p_ideal,p_ideal_filtered)

    plt.tick_params(axis='both', which='major', labelsize=int(_fontsize))
    plt.tick_params(axis='both', which='minor', labelsize=int(_fontsize))
    
    ax.set_xlim(0.8,p2[-1]+0.2)
    ax.set_ylim(0.8,ylim+0.2)
    
    plt.suptitle('',fontsize=12)
    
    fig = plt.gcf()
    matplotlib.pylab.legend(loc='best',fontsize='%d' % int(_fontsize))    
    DefaultSize = matplotlib.pylab.gcf().get_size_inches()
    fig.set_size_inches( (DefaultSize[0]/10, DefaultSize[1]/10) )
    fig.set_size_inches(7.25,7.25)
    
    if create_plot:
        plot_info = [''] * n_tables
        for i in range(0,n_tables):
            plot_info[i] = '+'.join(adapters[i])
            plot_info[i] = prefix[i] + '-' + plot_info[i]
            # make sure string is not too loong
            plot_info[i] = plot_info[i][:int(float(200)/float(n_tables))]
        plt.savefig('%s/%s.pdf' % ('.','_'.join(plot_info)), bbox_inches='tight')
        plt.savefig('%s/%s.png' % ('.','_'.join(plot_info)), bbox_inches='tight')
        print("PDF and PNG output written.")
    else:
        plt.show()
    return

########################################################################
# START OF THE PROGRAM
########################################################################
# For the meaning of the command line arguments, see the documentation
# of the function 'plot_multithreading_adapter_scaling' above.
help = '''
Creates a speedup plot based on the given tables containing adapter user and cpu times.
If multiple adapters are specified, then the cumulative user times are used to compute the speedup.
\n\n
Sample usage:\n
python plotmulticorespeedup.py -tables EulerFlow-p3' -prefixes \'EulerFlow-p3\' -legends \'2x Xeon  5-2650 @ 2.00GHz\' -ylim 16 -periteration -adapters \'Predictor+Corrector\' -cores 1 2 4 6 8 10 12 16 32 -hyperthreading
'''

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('-tables',nargs='+',required=True,help="A number of tables created via writecsvtable.py. The times read from the file in the first specified directory corresponding to the smallest thread count are considered as reference values.")
parser.add_argument('-prefixes',default=[''],nargs='+',required=True,help="Prefix of the files - usually the date of the test and an identifier for the machine and the MPI process that has written the output files. Must be supplied once per \'root_dir\' entry.")
parser.add_argument('-legends',nargs='+',required=True,help="Legend entry for the each data set - usually an identifier for the machine and the MPI process that has written the output files. Must be supplied once per \'root_dir\' entry.")
parser.add_argument('-adapters',default=['Total'],nargs='+',help="Name of the adapters separated by a \'+\'. (Use \'Total\' for the  cumulative time of all specified adapters. Does not make sense with \'per_iteration\' switched on.)")
parser.add_argument('-cores',nargs='+',required=True,help="The CPU cores (threads) you want to plot.")
parser.add_argument('-runs',default=1,help="Number of runs for all \'n\' and \'t\' combinations [default=1].")
parser.add_argument('-xticks',nargs='+',required=False,help="Ticks for the x-axis. Defaults to the thread counts 't'.")
parser.add_argument('-yticks',nargs='+',required=False,help="Ticks for the y-axis. Defaults to the x-axis ticks 'xticks'.")
parser.add_argument('-ylim',required=True,help="Upper limit for the y-axis.")
parser.add_argument('-hyperthreading', action='store_true', default=False,help="The last thread count corresponds to a hyperthreading run.")
parser.add_argument('-annotate', action='store_true', default=False,help="Annotate the plots with the speedup values.")
parser.add_argument('-periteration', action='store_true', default=False,help="Use the adapter times per iteration instead of the total times.")
parser.add_argument('-createplot', action='store_true', default=False,help="Creates a plot (PDF and PNG) in the working directory.")
parser.add_argument('-fontsize',default=10,required=False,help="Font size of the legend and tick labels. Axis labels are computed by ceiling the font size times a factor 1.2.")

args           = parser.parse_args();

tables         = args.tables
prefixes       = args.prefixes
legends        = args.legends
n_tables       = len(tables)
adapters       = [['']]*n_tables

for i in range(0,n_tables):
    adapters[i] = args.adapters[i].split('+')

cores          = args.cores
n_runs         = int(args.runs)
ylim           = float(args.ylim)
fontsize       = args.fontsize
hyperthreading = args.hyperthreading
annotate       = args.annotate
periteration   = args.periteration
createplot     = args.createplot

xticks         = args.xticks
if args.xticks is None:
    xticks     = cores
yticks         = args.yticks
if args.yticks is None:
    yticks     = xticks

plot_multithreading_adapter_scaling(tables,prefixes,legends,adapters,cores,n_runs,xticks,yticks,ylim,periteration,hyperthreading,annotate,createplot,fontsize)
