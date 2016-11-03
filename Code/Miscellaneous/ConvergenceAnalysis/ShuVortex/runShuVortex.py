#!/usr/bin/env python
#
# Run convergence tests.
#
#

import pandas as pd
from numpy import arange
import subprocess, os # batteries
shell = lambda shstring: subprocess.check_output(shstring, shell=True)

polyorders = arange(2,10)
width = 1. # make sure this is like in bash script.

depth = arange(1,6)
res = pd.DataFrame({'depth': depth})
res['meshsize'] = dx = width / 3.**depth  # actual meshsize we will get
res['numcells'] = width/dx
res['maxmeshsize'] = dx * 1.1 # as this is what the specfile wants

runner = "../RunScripts/runTemplatedSpecfile.sh"

settings = {}

settings['SIMBASE']="simulations/"
settings['ExaBinary']=shell("$(exa root)/$(exa find-binary EulerFlow)")
settings['ExaSpecfile']="ShuVortexConvergenceTpl.exahype"

# set initial data to use.
#settings['EXAHYPE_INITIALDATA']="MovingGauss2D"
settings['EXAHYPE_INITIALDATA']="ShuVortex"
# parameters for setting up the specfile
settings['EXASPEC_WIDTH']="15.0"
settings['EXASPEC_ENDTIME']="10.0"
# parameters deciding how frequently output is made. As a first criterion,
# 1 output dump with the highest resolution is 250MB.
settings['ExaConvOutputRepeat']="0.1"
settings['ExaVtkOutputRepeat']="0.5"

print "ExaHyPE convergence Analysis"
print "============================"

print "Can do convergence analysis on %fx%f sized domain with the following grid props: " % (width,width)
print res
print "Will do all these tests for these orders of the polynomial order:"
print polyorders

def start(polyorder, row):
	"""
	Starts a simulation

	@param polyorder is an integer
	@param row is a row in the res pdFrame
	"""
	env = os.environ.copy()
	env.update(settings)
	env['ExaMeshSize'] = str(row['maxmeshsize'])
	env["EXAREALMESHSIZE"] = str(row['meshsize'])
	env['ExapOrder'] = str(polyorder)
	env['ExaBinary'] = '{ExaBinary}-p{ExapOrder}'.format(**env)
	env['SIMDIR']="{SIMBASE}/p{ExapOrder}-meshsize{ExaMeshSize}/".format(**env)

	print "Starting p={ExapOrder}, maxmeshsize={ExaMeshSize}".format(**env)
	return subprocess.Popen([runner], env=env)

# adapt the number of runs / maximum cells for the polynomial degree
# adaptrunrange[<polyorder>] -> [<subset of maxmeshsizes>]
# numcells is array([   3.,    9.,   27.,   81.,  243.])
until = lambda maxcells: res[ res['numcells'] <= maxcells ]
adaptrunrange = { 2: until(243), 3: until(243), 4: until(243), 5: until(243),
		  6: until(81), 7: until(81),
		  8: until(27), 9: until(27) }

### TODO: Check command line arguments, idea:
###    -p  --polyorder => EXAPORDER
###    -m  --meshsize  => EXAMESHSIZE
### if not given, start all.
### Could also allow a range ==> See in my previous scripts
###                              about Neutron Star generation

for p, rows in adaptrunrange.iteritems():
	for i, row in rows.iterrows():
		start(p, row)

# script ends here.
# Watch the further output of the simulations on your computer using 'top' or 'htop'.
# Use 'multitail simulations/*/*.log' to watch what the simulations are doing.
# When finished (or even before), use 'finish-convergence-table.py' to inspect the results.

