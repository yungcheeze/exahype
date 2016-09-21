#!/usr/bin/env python
#
# Run convergence tests.
#
#

from numpy import arange
import subprocess # batteries

polyorder = 2. # ignored at this stage
depth = arange(2,6)
dx = (1./3.)**depth
numcells = 1./dx
reduced_maxmeshsizes = dx * 0.9

print "ExaHyPE convergence Analysis"
print "============================"

print "Will do convergence analysis with the following maximum mesh sizes: "
print dx
print "This is equivalent to these number of cells (unigrid, no AMR):"
print numcells

# run stuff in the background
procs = [
	subprocess.Popen(["./run-convergence.sh","-p",str(polyorder),"-m",str(maxmeshsize)])
	for maxmeshsize in reduced_maxmeshsizes
	]

# script ends here



