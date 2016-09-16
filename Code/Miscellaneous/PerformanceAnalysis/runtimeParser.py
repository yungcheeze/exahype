"""
.. module:: runtimeparser
  :platform: Unix, Windows, Mac
  :synopsis: Provides functions to read the CPU and user times of a set of Peano output files.
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>

:synopsis: Provides functions to read the CPU and user times of a set of Peano output files.
"""

import re


def max_level(filename):
  result = 0
  try:
    file = open(filename, "r" )
  except:
    print( "was not able to open file " + filename )
    return -1.0

  file.readline()
  for line in file:
    if re.search( "max-level=", line ):
      currentLevel = int( line.split( "max-level=" )[1].split(",")[0] )
      if currentLevel>result:
        result = currentLevel
  return result



def readColumnFromTable(filename,whichColumn):
  result = [] 
  try:
    file = open(filename, "r" )
  except:
    print( "was not able to open file " + filename )
    return []
  try:
    file.readline()
    for line in file:
      #if ( not re.search( "hline", line ) and ):
      #  print "got " + line
      try:
        data = float( line.split( "&" )[whichColumn] )
        result.append( data )
      except:
        result.append( 0.0 )
    return result
  except:
    print( "was not able to parse file " + filename + " holding a table ")
    return []


def getAdapterCountColumnFromTable(filename,adapterName):
  try:
    file = open(filename, "r" )
  except:
    print( "was not able to open file " + filename )
    return -1

  line = file.readline()
  counter = 0
  for col in line.split("&"):
    #print "study " + col.strip() + "  " 
    if col.strip()==adapterName.strip():
      return 2*counter-2
    counter = counter + 1
  return -1


def getAdapterRuntimeColumnFromTable(filename,adapterName):
  return getAdapterCountColumnFromTable(filename,adapterName) + 1


def parse_all_adapter_times(rootdir,prefix,process_counts,thread_counts,n_runs=1,cc='icpc',mode='TBB',per_iteration=False):
    """
    Reads multiple Peano output files with the naming pattern
       '<prefix>_n<#MPI proc.>_t<#threads/MPI proc.>_r<#run>_<cc>_<mode>.txt'

    and computes the average, minimum, and maximum CPU and user times for 
    each process and thread tuple. Averaging as well as searching for
    minima and maxima is performed over the number of runs 'n_runs'.

    Note:
       If the input array 'thread_counts' contains the entry 1,
       then the script looks for an output file
       with suffix 'None.txt' for this particular
       number of threads.

    Args:
       rootdir (str):
          Directory containing the Peano output files.
       prefix (str): 
          Prefix of the files - usually the date of the test and an identifier for the machine and the MPI process that has written the output files.
       process_counts (int[]):
          Number of MPI processes.
       thread_counts (int[]):
          Number of threads per node.
       n_runs (int):
          Number of runs for each 'n' and 'p' combination [default=1].
       cc (str):
          Compiler [default=gcc]. 
       mode (str):
          Shared memory mode [default=TBB].
       per_iteration (bool):
          Specifies if the times per iteration should be read in instead of the total times.

    Returns:
       A dict holding for each of the adapters a nested dict that holds the following key-value pairs:
          * 'n'           : (int)         Number of times this adapter was used.
          * 'avg_cputime' : (float[][]) CPU time spent within the adapter averaged over 'n_runs' runs.
          * 'min_cputime' : (float[][]) Minimum CPU time spent within the adapter out of 'n_runs' runs.
          * 'max_cputime' : (float[][]) Maximum CPU time spent within the adapter out of 'n_runs' runs.
          * 'avg_usertime': (float[][]) User time spent within the adapter averaged over 'n_runs' runs.
          * 'min_usertime': (float[][]) Minimum user time spent within the adapter out of 'n_runs' runs.
          * 'max_usertime': (float[][]) Maximum user time spent within the adapter out of 'n_runs' runs.
    """
    n_process_counts = len(process_counts)
    n_thread_counts  = len(thread_counts)
    
    result = { }
    for n in range(n_process_counts):
        for t in range(n_thread_counts):
            _mode = mode
            if int(thread_counts[t])==1:
                _mode = 'None'
            
            for r in range(n_runs):
                file_path = '%s/%s_n%s_t%s_r%d_%s_%s.txt' % (rootdir,prefix,str(process_counts[n]),str(thread_counts[t]),r+1,cc,_mode)
                result_r = parse_adapter_times(file_path,per_iteration)
                
                # We know all adapters after reading the first file and can thus initialise our data structure.
                if r==0 and t==0 and n==0:
                    for adapter in result_r:
                        result[adapter]                 = { }
                        result[adapter]['n']            = result_r[adapter]['n']
                        result[adapter]['avg_cputime']  = [[0.0]*n_thread_counts for i in range(n_process_counts)]
                        result[adapter]['min_cputime']  = [[float('inf')]*n_thread_counts for i in range(n_process_counts)]
                        result[adapter]['max_cputime']  = [[0.0]*n_thread_counts for i in range(n_process_counts)]
                        result[adapter]['avg_usertime'] = [[0.0]*n_thread_counts for i in range(n_process_counts)]
                        result[adapter]['min_usertime'] = [[float('inf')]*n_thread_counts for i in range(n_process_counts)]
                        result[adapter]['max_usertime'] = [[0.0]*n_thread_counts for i in range(n_process_counts)]
                
                for adapter in result:
                    result[adapter]['avg_cputime']  [n][t] += result_r[adapter]['cputime'] /float(n_runs);
                    result[adapter]['min_cputime']  [n][t]  = min(result[adapter]['min_cputime'] [n][t],result_r[adapter]['cputime'])
                    result[adapter]['max_cputime']  [n][t]  = max(result[adapter]['max_cputime'] [n][t],result_r[adapter]['cputime'])
                    result[adapter]['avg_usertime'] [n][t] += result_r[adapter]['usertime']/float(n_runs)
                    result[adapter]['min_usertime'] [n][t]  = min(result[adapter]['min_usertime'] [n][t],result_r[adapter]['usertime'])
                    result[adapter]['max_usertime'] [n][t]  = max(result[adapter]['max_usertime'] [n][t],result_r[adapter]['usertime'])
    
    return result

def parse_adapter_times(file_path,per_iteration=False):
    """
    Reads a single Peano output file and parses the user time spent within each adapter.
    
    Args:
       file_path (str):
          Path to the Peano output file.
       per_iteration (bool):
          Specifies if the times per iteration should be read in instead of the total times.

    Returns:
       A dict holding for each of the found adapters a nested dict that holds the following key-value pairs:
          * 'n'       : (int)    Number of times this adapter was used.
          * 'cputime' : (float) CPU time spent within the adapter.
          * 'usertime': (float) User time spent within the adapter.
    """
    result = { }
    try:
        file_handle=open(file_path)
        
        cputime_index  = 3
        usertime_index = 5
        if per_iteration:
            cputime_index  = 4
            usertime_index = 6
        
        for line in file_handle:
            # The trailing white space is important; do not remove it!
            anchor = '::repositories::RepositorySTDStack::logIterationStatistics() | '
            
            if anchor in line:
                segments = line.split('|')
                adapter = segments[1].strip();
                result[adapter]             = {}
                result[adapter]['n']        = int(segments[2].strip())
                result[adapter]['cputime']  = float(segments[cputime_index ].strip())
                result[adapter]['usertime'] = float(segments[usertime_index].strip())
    except:
        print ("Error: Could not process file '%s'!\n" % (file_path))
        raise
    return result

def sum_adapter_times(times,adapters,n_process_counts,n_thread_counts):
    """
    Sums the times for the specified adapters.
    
    Args:
       times ({}{}):
          A dict mapping to each of the adapters a nested dict that holds the following key-value pairs:
             * 'n'           : (int) Number of times this adapter was used.
             * 'avg_cputime' : (float[][]) CPU time spent within the adapter averaged over 'n_runs' runs.
             * 'min_cputime' : (float[][]) Minimum CPU time spent within the adapter out of 'n_runs' runs.
             * 'max_cputime' : (float[][]) Maximum CPU time spent within the adapter out of 'n_runs' runs.
             * 'avg_usertime': (float[][]) User time spent within the adapter averaged over 'n_runs' runs.
             * 'min_usertime': (float[][]) Minimum user time spent within the adapter out of 'n_runs' runs.
             * 'max_usertime': (float[][]) Maximum user time spent within the adapter out of 'n_runs' runs.
       adapters (str[]):
          The names of the adapters.

    Returns:
        A dict of the same structure as the nested dicts contained in
        argument 'times' holding cumulative times.
    """
    result                 = { }
    result['n']            = -1
    result['avg_cputime']  = [[0.0]*n_thread_counts for i in range(n_process_counts)]
    result['min_cputime']  = [[0.0]*n_thread_counts for i in range(n_process_counts)]
    result['max_cputime']  = [[0.0]*n_thread_counts for i in range(n_process_counts)]
    result['avg_usertime'] = [[0.0]*n_thread_counts for i in range(n_process_counts)]
    result['min_usertime'] = [[0.0]*n_thread_counts for i in range(n_process_counts)]
    result['max_usertime'] = [[0.0]*n_thread_counts for i in range(n_process_counts)]

    for n in range(n_process_counts):
        for t in range(n_thread_counts):
            for adapter in times:
                if adapter in adapters:
                    result['n']                    = times[adapter]['n']
                    result['avg_cputime']  [n][t] += times[adapter]['avg_cputime'] [n][t]
                    result['min_cputime']  [n][t] += times[adapter]['min_cputime'] [n][t]
                    result['max_cputime']  [n][t] += times[adapter]['max_cputime'] [n][t]
                    result['avg_usertime'] [n][t] += times[adapter]['avg_usertime'][n][t]
                    result['min_usertime'] [n][t] += times[adapter]['min_usertime'][n][t]
                    result['max_usertime'] [n][t] += times[adapter]['max_usertime'][n][t]
    
    return result

def sum_all_adapter_times(times,n_process_counts,n_thread_counts):
    """
    Computes the sum of the CPU and user times over all adapters.
    
    Args:
       times ({}{}):
          A dict mapping to each of the adapters a nested dict that holds the following key-value pairs:
             * 'n'           : (int)       Number of times this adapter was used.
             * 'avg_cputime' : (float[][]) CPU time spent within the adapter averaged over 'n_runs' runs.
             * 'min_cputime' : (float[][]) Minimum CPU time spent within the adapter out of 'n_runs' runs.
             * 'max_cputime' : (float[][]) Maximum CPU time spent within the adapter out of 'n_runs' runs.
             * 'avg_usertime': (float[][]) User time spent within the adapter averaged over 'n_runs' runs.
             * 'min_usertime': (float[][]) Minimum user time spent within the adapter out of 'n_runs' runs.
             * 'max_usertime': (float[][]) Maximum user time spent within the adapter out of 'n_runs' runs.

    Returns:
       A dict holding the following key-value pairs:
          * 'avg_cputime' : (float[][]) Total CPU time averaged over 'n_runs' runs.
          * 'min_cputime' : (float[][]) Total minimum CPU time out of 'n_runs' runs.
          * 'max_cputime' : (float[][]) Total maximum CPU time out of 'n_runs' runs.
          * 'avg_usertime': (float[][]) Total user time averaged over 'n_runs' runs.
          * 'min_usertime': (float[][]) Total minimum user time spent out of 'n_runs' runs.
          * 'max_usertime': (float[][]) Total maximum user time spent out of 'n_runs' runs.
    """
    result = { }
    result['n']            = -1
    result['avg_cputime']  = [[0.0]*n_thread_counts for i in range(n_process_counts)]
    result['min_cputime']  = [[0.0]*n_thread_counts for i in range(n_process_counts)]
    result['max_cputime']  = [[0.0]*n_thread_counts for i in range(n_process_counts)]
    result['avg_usertime'] = [[0.0]*n_thread_counts for i in range(n_process_counts)]
    result['min_usertime'] = [[0.0]*n_thread_counts for i in range(n_process_counts)]
    result['max_usertime'] = [[0.0]*n_thread_counts for i in range(n_process_counts)]

    for n in range(n_process_counts):
        for t in range(n_thread_counts):
            for adapter in times:
                result['avg_cputime']  [n][t] += times[adapter]['avg_cputime'] [n][t]
                result['min_cputime']  [n][t] += times[adapter]['min_cputime'] [n][t]
                result['max_cputime']  [n][t] += times[adapter]['max_cputime'] [n][t]
                result['avg_usertime'] [n][t] += times[adapter]['avg_usertime'][n][t]
                result['min_usertime'] [n][t] += times[adapter]['min_usertime'][n][t]
                result['max_usertime'] [n][t] += times[adapter]['max_usertime'][n][t]
    
    return result
