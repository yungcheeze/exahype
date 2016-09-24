#!/usr/bin/env python
#
# Run convergence tests.
#
#

from numpy import arange
import subprocess # batteries

polyorder = 3. # make sure binary is compiled with this order
width = 15. # make sure this is like in bash script.

depth = arange(2,6)
dx = width / 3.**depth
numcells = width/dx
reduced_maxmeshsizes = dx * 1.1

def listprint(l):
	for t in enumerate(l): print "Run %d: %f"%t

print "ExaHyPE convergence Analysis"
print "============================"

print "Will do convergence analysis on %fx%f sized domain with the following maximum mesh sizes: " % (width,width)
listprint(dx)
print "With the reduced mesh sizes"
listprint(reduced_maxmeshsizes)
print "This is equivalent to these number of cells (unigrid, no AMR):"
listprint(numcells)

#import sys; sys.exit(0)

# run stuff in the background
procs = [
	subprocess.Popen(["./run-convergence.sh","-p",str(polyorder),"-m",str(maxmeshsize)])
	for maxmeshsize in reduced_maxmeshsizes
	]

# script ends here.
# Watch the further output of the simulations on your computer using 'top' or 'htop'.
# Use 'multitail simulations/*/*.log' to watch what the simulations are doing.
# When finished (or even before), use 'finish-convergence-table.py' to inspect the results.
