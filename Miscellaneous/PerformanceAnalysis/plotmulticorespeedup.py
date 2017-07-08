#!/usr/bin/python3
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
from argparse import RawTextHelpFormatter

import matplotlib.pylab
import matplotlib.pyplot as plt

import operator

import csv

def compute_speedup(t,t_ref):
    """
    Args:
       t  (float[]):
          Array holding times.
       t_ref (float):
          Reference time.
    Returns:
       Speedup (float[]) with respect to 't_ref'.
    """
    n_t = len(t)
    S    = [0.0]*n_t
    
    for i in range(0,n_t):
        S[i] = t_ref/t[i]
    return S

def plot_scaling(ax,x,y,label_,color_="blue",marker_="s",markerfacecolor_="blue",hyperthreading=False,annotate=False):
    """
    Plots y over x. If 'hyperthreading' is selected, the last value is plotted 
    separately. If 'annotate' is selected, the plot is annotated with the 'y'
    values.
    """
    x2 = x
    
    # In the hyperthreading case, we only want to plot a line
    # for the real threads. The hyperthreading value is plotted separately
    if not hyperthreading:
        plt.plot(x,y,label=r"%s" % label_,markersize=8,marker=marker_,markevery=1,lw=1.1,color=color_,markerfacecolor=markerfacecolor_,markeredgecolor=color_) 
    else:
        x2 = x[0:-1] # [inclusive:exclusive]
        x2.append(x2[-1]+2)
        
        plt.plot(x2[0:-1],y[0:-1],label=r"%s" % label_,markersize=8,marker=marker_,markevery=1,lw=1.1,color=color_,markerfacecolor=markerfacecolor_,markeredgecolor=color_) 
        plt.plot(x2[-1],y[-1],label="",markersize=8,marker=marker_,markevery=1,lw=1.1,color=color_,markerfacecolor=markerfacecolor_,markeredgecolor=color_) 
    
    if annotate:
        for i,j in zip(x2,y):
            ax.annotate("%.2f" % j,xy=(i,j+0.5),color="blue")
    return


def query_table(filename,mesh,order,adapter,cores): 
    '''
    Read a table and filter out certain
    columns 
    Args:
    filename
       name of the csv file
    order
       approximaation order
    adapter
       the adapter
    cores
       list of cores
    '''
    print(filename)
    print(mesh)
    print(order)
    datafile    = open(filename, 'r')
    reader      = csv.reader(datafile,delimiter='&')
    sorted_data = list(filter(lambda x : x[0]==mesh and x[1]==order and x[5]==adapter and x[8] in cores, reader))
    datafile.close() 
    return sorted_data


'''
.. module:: usertimeplot
  :platform: Unix, Windows, Mac
  :synopsis: Creates a speedup plot based on Peano output files with specific file naming pattern.
   
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>

:synopsis: Creates a speedup plot based on user and CPU times stored in CSV tables.
'''

def plot_multithreading_adapter_scaling(table,prefix,mesh,legend,order,adapters,cores,xticks,yticks,ylim,per_iteration=False,hyperthreading=False,annotate=False,create_plot=False,_fontsize=10):
    '''
    Creates a scaling plot for the cumulative user time spent within the 
    specified adapters.
   
    Args:
      table (str):
         A csv file
      prefix (str):
         prefix for the output file
      mesh (str):
        the mesh
      legend (str[]):
         Legend entry for the each data set - usually an identifier for the machine and the MPI process that has written the output files. Must be supplied once per 'root_dir' entry.
      order (str[]):
         Approximation orders
      adapter (str[]):
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
    
    n_legends = len(legend)
    n_orders  = len(order)
    n_cores   = len(cores)
    
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
    for i in range(0,n_legends):
        for N in order:
            cumulative_user_times = [0.0]*n_cores
            cumulative_cpu_times  = [0.0]*n_cores
            iterations = 1.0
            for adapter in adapters[i]:
                data = query_table(table,mesh,N,adapter,cores)
                #for row in data:
                #    print(row)
                cumulative_user_times = list(map(lambda x,y : x+y, cumulative_user_times, list(map(float,[row[11] for row in data]))))
                cumulative_cpu_times  = list(map(lambda x,y : x+y, cumulative_cpu_times,  list(map(float,[row[12] for row in data]))))
                
                if periteration:
                    iterations               = float(data[0][10]) # must be the same for all cores
                    cumulative_user_times[:] = [x / iterations for x in cumulative_user_times]
                    cumulative_cpu_times[:]  = [x / iterations for x in cumulative_cpu_times]
                    
            print(cumulative_user_times)
            print(cumulative_cpu_times)

            if i==0:
                t_ref = cumulative_user_times[0]
            
            speedup_measured = compute_speedup(cumulative_user_times,t_ref)

            # Measured speedup
            plot_scaling(ax,p,speedup_measured,legend[i]+"(N="+N+")",colors[i],markers[i],markerfacecolors[i],hyperthreading,annotate)

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
        plot_info = [''] * n_legends
        for i in range(0,n_legends):
            plot_info[i] = '+'.join(adapters[i])
            # make sure string is not too loong
            plot_info[i] = plot_info[i][:int(float(200)/float(n_legends))]
        plt.savefig('./%s-%s.pdf' % (prefix,'_'.join(plot_info)), bbox_inches='tight')
        plt.savefig('./%s-%s.png' % (prefix,'_'.join(plot_info)), bbox_inches='tight')
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
parser.add_argument('-table',required=True,help="A csv file with & separator.")
parser.add_argument('-prefix',required=True,help="Prefix of the plot files.")
parser.add_argument('-mesh',required=True,help="Mesh identifier.")
parser.add_argument('-legend',nargs='+',required=True,help="Legend entry for the each data set - usually an identifier for the machine and the MPI process that has written the output files. Must be supplied once per \'root_dir\' entry.")
parser.add_argument('-adapter',default=['Total'],nargs='+',help="Name of the adapters separated by a \'+\'. (Use \'Total\' for the  cumulative time of all specified adapters. Does not make sense with \'per_iteration\' switched on.)")
parser.add_argument('-order',nargs='+',required=True,help="The orders/patch sizes you want to plot")
parser.add_argument('-cores',nargs='+',required=True,help="The CPU cores (threads) you want to plot.")
parser.add_argument('-xticks',nargs='+',required=False,help="Ticks for the x-axis. Defaults to the thread counts 't'.")
parser.add_argument('-yticks',nargs='+',required=False,help="Ticks for the y-axis. Defaults to the x-axis ticks 'xticks'.")
parser.add_argument('-ylim',required=True,help="Upper limit for the y-axis.")
parser.add_argument('-hyperthreading', action='store_true', default=False,help="The last thread count corresponds to a hyperthreading run.")
parser.add_argument('-annotate', action='store_true', default=False,help="Annotate the plots with the speedup values.")
parser.add_argument('-periteration', action='store_true', default=False,help="Use the adapter times per iteration instead of the total times.")
parser.add_argument('-createplot', action='store_true', default=False,help="Creates a plot (PDF and PNG) in the working directory.")
parser.add_argument('-fontsize',default=10,required=False,help="Font size of the legend and tick labels. Axis labels are computed by ceiling the font size times a factor 1.2.")

args           = parser.parse_args();

table          = args.table
prefix         = args.prefix
mesh           = args.mesh
legend         = args.legend
n_legends      = len(legend)
adapters       = [['']]*n_legends

for i in range(0,n_legends):
    adapters[i] = args.adapter[i].split('+')

order          = args.order
cores          = args.cores
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
    
print(yticks)

plot_multithreading_adapter_scaling(table,prefix,mesh,legend,order,adapters,cores,xticks,yticks,ylim,periteration,hyperthreading,annotate,createplot,fontsize)
