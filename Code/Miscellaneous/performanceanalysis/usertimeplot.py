import argparse
from argparse import RawTextHelpFormatter

import matplotlib.pylab
import matplotlib.pyplot as plt

import runtimeparser as rp
import hpclib as hpc
from plotting import scalingplot as sp

'''
.. module:: usertimeplot
  :platform: Unix, Windows, Mac
  :synopsis: Creates a speedup plot based on Peano output files with specific file naming pattern.
   
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>

:synopsis: Creates a speedup plot based on  Peano output files with specific file naming pattern.
'''

def plot_multithreading_adapter_scaling(root_dir,prefix,legend,adapters,process_counts,thread_counts,n_runs,cc,mode,xticks,yticks,ylim,per_iteration=False,hyperthreading=False,annotate=False,create_pdf=False,_fontsize=10):
    '''
    Creates a scaling plot for the cumulative user time spent within the 
    specified adapters.
   
    Args:
      root_dir (str[]):
         Directories containing the Peano output files. (Implementation does currently only support one root dir element.)
      prefix (str[]):
         Prefix of the files - usually the date of the test and an identifier for the machine and the MPI process that has written the output files. Must be supplied once per 'root_dir' entry.
      legend (str[]):
         Legend entry for the each data set - usually an identifier for the machine and the MPI process that has written the output files. Must be supplied once per 'root_dir' entry.
      adapters (str[]):
         Name of the adapters. (Use 'Total' for the  cumulative time forall adapters. Does not make sense with per_iteration switched on.)
      process_counts (int[]):
         MPI process counts.
      thread_counts  (int[]):
         Threads per MPI process.
      n_runs (int):
         Number of runs for each 'n' and 't' combination [default=1].
      cc (str): 
         Compiler.
      mode (str):
         Shared memory mode.
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
    n_root_dir       = len(root_dir)
    
    colors           = ['k','k','k','k','k','k']
    markerfacecolors = ['None','None','None','k','k','k']
    markers          = ['o','^','s','o','^','s']
    
    n_process_counts = len(process_counts)
    n_thread_counts  = len(thread_counts)
    
    # Plotting
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    
    p       = map(int,thread_counts)
    p2      = p
    p_ticks = map(str,p)
    
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
    
    for i in range(0,n_root_dir):
        # Read in user and CPU times
        times            = rp.parse_all_adapter_times(root_dir[i],prefix[i],process_counts,thread_counts,n_runs,cc[0],mode[0],per_iteration)
        total_times      = rp.sum_all_adapter_times(times,n_process_counts,n_thread_counts)
        times['Total']   = total_times
        cumulative_times = rp.sum_adapter_times(times,adapters[i],n_process_counts,n_thread_counts)
        
        if i==0:
            t_ref = cumulative_times['avg_usertime'][0][0]
        
        speedup_measured = hpc.compute_speedup_2(cumulative_times['avg_usertime'],t_ref)
        speedup_measured = speedup_measured[0];
    
        # Measured speedup
        sp.plot_scaling(ax,p,speedup_measured,legend[i],colors[i],markers[i],markerfacecolors[i],hyperthreading,annotate)

    plt.ylabel(r'speedup',         fontsize=float(1.2*float(_fontsize)))    
    plt.xlabel(r'number of cores', fontsize=float(1.2*float(_fontsize)))
    plt.grid(True)

    p_ticks_filtered = p_ticks
    for ip in range(0,len(p_ticks)):
        if 'HT' not in p_ticks_filtered[ip]:
            if p_ticks_filtered[ip] not in xticks:
                p_ticks_filtered[ip] = ''
    plt.xticks(p2,p_ticks_filtered)

    p_ideal_filtered = map(int,p_ideal)
    p_ideal_filtered = map(str,p_ideal_filtered)
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
    
    if create_pdf:
        plot_info = [''] * n_root_dir
        for i in range(0,n_root_dir):
            plot_info[i] = '+'.join(adapters[i])
            plot_info[i] = prefix[i] + '-' + plot_info[i]
        plt.savefig('%s/%s.pdf' % ('.','_'.join(plot_info)), bbox_inches='tight')
        print("PDF output written.")
    plt.show()
    return

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
parser.add_argument('-path',nargs='+',required=True,help="Directories containing the Peano output files. The times read from the file in the first specified directory corresponding to the smallest thread count are considered as reference values.")
parser.add_argument('-prefix',default=[''],nargs='+',required=True,help="Prefix of the files - usually the date of the test and an identifier for the machine and the MPI process that has written the output files. Must be supplied once per \'root_dir\' entry.")
parser.add_argument('-legend',nargs='+',required=True,help="Legend entry for the each data set - usually an identifier for the machine and the MPI process that has written the output files. Must be supplied once per \'root_dir\' entry.")
parser.add_argument('-adapter',default=['Total'],nargs='+',help="Name of the adapters separated by a \'+\'. (Use \'Total\' for the  cumulative time of all specified adapters. Does not make sense with \'per_iteration\' switched on.)")
parser.add_argument('-n',default=[1],nargs='+',help="MPI process counts [default=1].")
parser.add_argument('-t',default=[1],nargs='+',required=True,help="Threads per MPI process [default=1].")
parser.add_argument('-r',default=1,help="Number of runs for all \'n\' and \'t\' combinations [default=1].")
parser.add_argument('-cc',default='icpc',nargs='+',help="Compiler [default=\'icpc\']")
parser.add_argument('-mode',default='TBB',nargs='+',help="Shared memory mode [default=\'TBB\']")
parser.add_argument('-xticks',nargs='+',required=False,help="Ticks for the x-axis. Defaults to the thread counts 't'.")
parser.add_argument('-yticks',nargs='+',required=False,help="Ticks for the y-axis. Defaults to the x-axis ticks 'xticks'.")
parser.add_argument('-ylim',required=True,help="Upper limit for the y-axis.")
parser.add_argument('-hyperthreading', action='store_true', default=False,help="The last thread count corresponds to a hyperthreading run.")
parser.add_argument('-annotate', action='store_true', default=False,help="Annotate the plots with the speedup values.")
parser.add_argument('-per_iteration', action='store_true', default=False,help="Use the adapter times per iteration instead of the total times.")
parser.add_argument('-create_pdf', action='store_true', default=False,help="Creates a PDF plot in the working directory.")
parser.add_argument('-fontsize',default=10,required=False,help="Font size of the legend and tick labels. Axis labels are computed by ceiling the font size times a factor 1.2.")

args           = parser.parse_args();

root_dir       = args.path
prefix         = args.prefix
legend         = args.legend
n_root_dir     = len(root_dir)
adapter        = [['']]*n_root_dir

for i in range(0,n_root_dir):
    adapter[i] = args.adapter[i].split('+')

process_counts = args.n
thread_counts  = args.t
n_runs         = int(args.r)
cc             = args.cc
mode           = args.mode
ylim           = float(args.ylim)
fontsize       = args.fontsize
hyperthreading = args.hyperthreading
annotate       = args.annotate
per_iteration  = args.per_iteration
create_pdf     = args.create_pdf

xticks         = args.xticks
if args.xticks is None:
    xticks     =thread_counts
yticks         = args.yticks
if args.yticks is None:
    yticks     =xticks

plot_multithreading_adapter_scaling(root_dir,prefix,legend,adapter,process_counts,thread_counts,n_runs,cc,mode,xticks,yticks,ylim,per_iteration,hyperthreading,annotate,create_pdf,fontsize)
