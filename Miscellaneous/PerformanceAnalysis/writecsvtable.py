#!/user/bin/python3
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
import argparse
from argparse import RawTextHelpFormatter

import re
import os
import csv

import operator

import runtimeParser

'''
.. module:: extractCSVTable
  :platform: Unix, Windows, Mac
  :synopsis: Extracts performance metrics from Peano output files with specific file naming pattern.
   
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>

:synopsis: Extracts performance metrics from Peano output files with specific file naming pattern.
'''

def extract_table(root_dir,prefix):
    '''
    Extracts performance metrics from Peano output files with specific file naming pattern.
   
    Args:
      root_dir (str):
         Directory containing the Peano output files.
      prefix (str):
         Prefix of the files - usually the date of the test and an identifier for the test.
    '''
 
    # collect filenames
    with open(prefix+'.csv', 'w') as csvfile:
        csvwriter = csv.writer(csvfile, delimiter='&',quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for filename in os.listdir(root_dir):
            if filename.endswith(".out") and filename.startswith(prefix):
                match = re.search(prefix+'-n([0-9]+)-t([0-9]+)-c([0-9]+)-([A-Za-z]+)-([A-Za-z]+)\.out',filename)
                nodes = match.group(1)
                tasks = match.group(2)
                cores = match.group(3)
                mode  = match.group(4)
                cc    = match.group(5)
                    
                times = runtimeParser.parse_adapter_times(filename) 
                
                for adapter in times:
                    iterations = times[adapter]['n']
                    usertime   = times[adapter]['usertime']
                    cputime    = times[adapter]['cputime']
 
                    csvwriter.writerow([nodes,tasks,cores,adapter,iterations,usertime,cputime,cc,mode])

def sort_table(filename):
    '''
    Sorts the rows of the file according to nodes,tasks,cores,adapter name.
    See: https://stackoverflow.com/a/17109098
    '''
    datafile    = open(filename, 'r')
    reader      = csv.reader(datafile,delimiter='&')
    sorted_data = sorted(reader, key=lambda x: (x[3],int(x[0]),int(x[1]),int(x[2])))
    datafile.close() 
 
    with open(filename, 'w') as datafile:
        writer = csv.writer(datafile, delimiter='&',quotechar='|', quoting=csv.QUOTE_MINIMAL)
        writer.writerows(sorted_data)

########################################################################
# START OF THE PROGRAM
########################################################################
# For the meaning of the command line arguments, see the documentation
# of the function 'extrat_table' above.
help = '''
Extract performance metrics from Peano output files with specific file naming pattern
and write them to a csv file with name
<prefix_<cc>_<mode>.csv

\n\n
Sample usage:\n
python extract-table-dominic.py -path \'examples/151217_phi1_node/'
'''

parser = argparse.ArgumentParser(description=help,formatter_class=RawTextHelpFormatter)
parser.add_argument('-path',required=True,help="Directory containing the Peano output files.")
parser.add_argument('-prefix',required=True,help="Prefix of the Peano output files.")

args     = parser.parse_args();

root_dir = args.path
prefix   = args.prefix

extract_table(root_dir,prefix)
sort_table(root_dir+"/"+prefix+".csv")
