#!/usr/bin/env python
#
# Try this script with
#   ipython -i finish-convergence-table.py
# to inspect all the lists and dicts.
#
# -- SK, 2016

from numpy import genfromtxt
# batteries:
from glob import glob
from os import path, stat
from re import match

# which quantity shall we look at, in the moment?
quantityfile="output/error-rho.asc"
simulations = glob('simulations/*')

def read_simulation_params(envfile):
	"""
	Small function to parse files in this format:
	EXAMESHSIZE=0.0037037037037
	EXAPORDER=2.0
	EXAHYPE_INITIALDATA=MovingGauss2D
	...
	Returns dictionary like
	{
	'EXAMESHSIZE': '0.0037037037037',
	'EXAPORDER': '2.0',
	...
	}
	"""
	params = {}
	with open(envfile, 'r') as fh:
		for line in fh:
			res = match(r'^([a-zA-Z0-9]+)=(.*)$', line)
			if res:
				params[ res.group(1) ] = res.group(2)
	return params

# just a helper function to detect empty files
filenotempty = lambda f: stat(f).st_size != 0
simfile = lambda fname: lambda simdir: path.join(simdir, fname)

# fault tolerant genfromtxt
def trygenfromtxt(f, *args, **kwargs):
	try:
		return genfromtxt(f, *args, **kwargs)
	except e:
		return None

# filter out empty simulations
simulations = filter(filenotempty, simulations)

print "Have read %d non-empty simulations: " % len(simulations), simulations

# determine filenames
quantityfiles = map(simfile(quantityfile), simulations)
paramfiles = map(simfile('parameters.env'), simulations)

# read files
params = map(read_simulation_params, paramfiles)

# genfromtxt with header detection is very sensible to the first line's format.
# this will not work:
#    # 1:plotindex ,2:time ,3:l1norm ,4:l2norm ,5:max ,6:min ,7:avg ,
# in constrast, the line has to look like
#    1:plotindex 2:time 3:l1norm 4:l2norm 5:max 6:min 7:avg
# or even better
#    plotindex time l1norm l2norm max min avg
errors = [trygenfromtxt(f, names=True) for f in quantityfiles]

# now all data is available:
#
# params[0]['EXAMESHSIZE']
# errors[0]['l1norm']
# etc.


