import argparse
from argparse import RawTextHelpFormatter

import re

import runtimeParser as rp
import hpclib as hpc

'''
.. module:: usertimeplot
  :platform: Unix, Windows, Mac
  :synopsis: Extracts performance metrics from Peano output files with specific file naming pattern.
   
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>

:synopsis: Extracts performance metrics from Peano output files with specific file naming pattern.
'''

def extract_table(root_dir,prefix,nodes,tasks,cores,runs,cc,mode)
    '''
    Extracts performance metrics from Peano output files with specific file naming pattern.
   
    Args:
      root_dir (str):
         Directory containing the Peano output files.
      prefix (str):
         Prefix of the files - usually the date of the test and an identifier for the test.
    '''    
   	# collect filenames
   	for filename in os.listdir(root_dir):
    if filename.endswith(".out") and filename.startswith(prefix) 
        match = re.search(prefix+'-n([0-9]+)-t([0-9]+)-c([0-9]+)-([A-Za-z]+)-([A-Za-z]+)\.out')
        nodes = match.group(1)
        tasks = match.group(2)
        cores = match.group(3)
        mode  = match.group(4)
        cc    = match.group(5)
        
        print("Found file "+filename)
        print("Extracted the following metadata from file name:")
        print("nodes="nodes)
        print("tasks="tasks)
        print("cores="cores)
        print("mode="mode)
        print("cc="cc)
        continue
    else:
        print("Did not find any file with prefix "+prefix+" and suffix \'.out\'.")

########################################################################
# START OF THE PROGRAM
########################################################################
# For the meaning of the command line arguments, see the documentation
# of the function 'extrat_table' above.
help = '''
Extra performance metrics from Peano output files with specific file naming pattern
and write them to a csv file with name
<cc>_<TBB>

\n\n
Sample usage:\n
python extract-table-dominic.py -path \'examples/151217_phi1_node/'
'''

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('-path',required=True,help="Directories containing the Peano output files. The times read from the file in the first specified directory corresponding to the smallest thread count are considered as reference values.")

args           = parser.parse_args();

root_dir       = args.path

extract_table(root_dir)
