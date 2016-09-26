#!/usr/bin/env python
#
# Run convergence tests.
#
#

import pandas as pd
from numpy import arange
import subprocess, os # batteries

polyorders = arange(2,10)
width = 15. # make sure this is like in bash script.

depth = arange(1,6)
res = pd.DataFrame({'depth': depth})
res['meshsize'] = dx = width / 3.**depth  # actual meshsize we will get
res['numcells'] = width/dx
res['maxmeshsize'] = dx * 1.1 # as this is what the specfile wants

print "ExaHyPE convergence Analysis"
print "============================"

print "Can do convergence analysis on %fx%f sized domain with the following grid props: " % (width,width)
print res
print "Will do all these tests for these orders of the polynomial order:"
print polyorders

def start(polyorder, row):
	# row in res
	# somehow the calling convention for run-convergence.sh is shitty.
	env = os.environ.copy()
	env["EXAREALMESHSIZE"] = str(row['meshsize'])
	print "Starting p=%d, maxmeshsize=%f" % (polyorder, row['maxmeshsize'])
	return subprocess.Popen(["./run-convergence.sh","-p",str(polyorder),"-m",str(row['maxmeshsize'])], env=env)

# adapt the number of runs / maximum cells for the polynomial degree
# adaptrunrange[<polyorder>] -> [<subset of maxmeshsizes>]
# numcells is array([   3.,    9.,   27.,   81.,  243.])
until = lambda maxcells: res[ res['numcells'] <= maxcells ]
adaptrunrange = { 2: until(243), 3: until(243), 4: until(243), 5: until(243),
		  6: until(81), 7: until(81),
		  8: until(27), 9: until(27) }

for p, rows in adaptrunrange.iteritems():
	for i, row in rows.iterrows():
		start(p, row)

# script ends here.
# Watch the further output of the simulations on your computer using 'top' or 'htop'.
# Use 'multitail simulations/*/*.log' to watch what the simulations are doing.
# When finished (or even before), use 'finish-convergence-table.py' to inspect the results.
