#!/usr/bin/env python
#
# Try this script with
#   ipython -i finish-convergence-table.py
# to inspect all the lists and dicts.
#
# -- SK, 2016

import numpy as np
import pandas as pd
# batteries:
from glob import glob
from os import path, stat
from re import match, sub
from os import path
from string import Template # HTML output
# just for statistics
from datetime import datetime
from socket import gethostname
from getpass import getuser

# which quantity shall we look at, in the moment?
quantity = 'error-rho.asc'
quantityfile="output/" + quantity
simulations = glob('simulations/p3*/') # mind the trailing slash to glob only directories
report_templatefile = "report/report-template.html"
report_outputfile = "simulations/generated-report.html"

def read_simulation_params(envfile):
	"""
	Small function to parse files in this format:
	EXAMESHSIZE=0.0037037037037
	EXAPORDER=2.0
	EXAHYPE_INITIALDATA=MovingGauss2D
	...
	Returns dictionary like
	{
	'EXAMESHSIZE': 0.0037037037037,
	'EXAPORDER': 2.0,
	'EXAHYPE_INITIALDATA': 'MovingGauss2D',
	...
	}
	"""
	params = {}
	with open(envfile, 'r') as fh:
		for line in fh:
			res = match(r'^([a-zA-Z0-9_]+)=(.*)$', line)
			if res:
				k, v = res.group(1), res.group(2)
				# try to cast the value as a float,
				# as this is a common data type in our buisness
				try:
					params[k] = float(v)
				except:
					params[k] = v
	return params

# just a helper function to detect empty files
filenotempty = lambda f: stat(f).st_size != 0
simfile = lambda fname: lambda simdir: path.join(simdir, fname)

# fault tolerant genfromtxt
def trygenfromtxt(f, *args, **kwargs):
	try:
		return np.genfromtxt(f, *args, **kwargs)
	except:
		return None

# filter out empty simulations
# @todo, does not work right now, have to check on `quantityfiles`
simulations = filter(filenotempty, simulations)

print "Have read %d non-empty simulations: " % len(simulations)
for t in enumerate(simulations): print " %i. %s" % t

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

## Numpy version:
# errors = [trygenfromtxt(f, names=True) for f in quantityfiles]

## Pandas version: a list of data frames
errors = [pd.read_csv(qf, delim_whitespace=True) for qf in quantityfiles]

# get the meshsizes of the individual simulations
meshsizes = np.array([ pi['EXAMESHSIZE'] for pi in params])
# something like [90, 10, 30]
ncells = 1. / meshsizes

## we can also convert the parameters to a pandas table
paramtable = pd.DataFrame(params)
# beautify the parameter table: keep only columns with "EXA" inside
cols = [col for col in paramtable.columns if 'EXA' in col]
paramtable = paramtable[cols]
# beautify the parameter table: shorten the paths
possiblepathcolumns = ['EXABINARY', 'EXASPECFILE']
for col in possiblepathcolumns:
	try: paramtable[col] = paramtable[col].map(path.basename)
	except e: pass


# merge data frames to multiindex
idxcells = 'nCells' # name of 'number of cells' column
errors = pd.concat(errors, keys=ncells, names=[idxcells])

# how much data rows can we compare? This is the number of maximum
# common evolution steps
## Numpy version:
#maxcomparisons = min([len(dataset) for dataset in errors])
## Pandas version:
maxcomparisons = errors.groupby(level=idxcells).size().min()

print "Have this error table from all simulations, with the maximum"
print "comparable number of rows = %d" % maxcomparisons
print errors

# either strip down datasets to common errors and analyse all of them
## Numpy version:
#ce = [dataset[:maxcomparisons] for dataset in errors]

# or, for simplicity for the time being, choose the LAST time evolution step
# as the most meaningful in terms of accuracy.
# This means `celast` is an array of single row structured numpy arrays.
## Numpy version:
#celast = [ dataset[maxcomparisons-1] for dataset in errors ]
## Pandas version:

idxplotindex = 'plotindex' # the column counting the rows in each simulation
ceslice = errors.query('%s <= %d' % (idxplotindex, maxcomparisons))
celast = ceslice.groupby(level=idxcells).tail(1)

# common errors, the dataset we will work with now.
ce = ceslice

# to get rid of the index:
ce = ce.reset_index()
ce.sort(idxcells, inplace=True) # inserts a "level_1" column

# find the reference nCells against which we make the convergence test.
# This is usually the one nCell at the same plotindex smaller than the
# actual one.
def findSmallerNcellsThan(num):
	smallerNcells = ncells[ ncells < num ]
	return None if not len(smallerNcells) else smallerNcells.max()

idxprev = 'nSmaller' # column name of linked idxcells ('nCells')
ce[idxprev] = ce[idxcells].apply(findSmallerNcellsThan)

# insert the convergence measures for these columns
columns = 'l1norm l2norm max min avg'.split()
outcols = ["o"+c for c in columns]
idxtime = 'time' # the column for time

# add the new columns 
for newcol in outcols:
	ce[newcol] = np.zeros(len(ce))

# compute the actual convergence number for each column
for rowindex, row in ce.iterrows():
	print "I am in row %d and row[%s] is %s" % (rowindex, idxprev, str(row[idxprev]))

	if not row[idxprev] or np.isnan(row[idxprev]):
		# there is no smaller resolution available
		continue

	targetrow = ce[(ce[idxcells] == row[idxprev]) & (ce[idxplotindex] == row[idxplotindex])]
	if targetrow.empty:
		# could not find a previous step, so silently ignore
		print "Did not find a target row"
		continue

	if len(targetrow) != 1:
		print "Could not find unique counterpart for row ", row
		print "Found instead %d rows: " % len(targetrow), targetrow
		continue

	for col, outcol in zip(columns, outcols):
		value = np.log(row[col] / targetrow[col]) / np.log( targetrow[idxcells] / row[idxcells] )
		print "Computed row[%s]=%f" % (outcol, value)
		ce.set_value(rowindex, outcol, value)

# print out a subset of the table
print "Comptuted this convergence table for the individual reductions"
print "(as l1norm, infnorm=max, etc.)"
convergence_table = ce[[idxcells,idxprev,idxplotindex] + outcols]
print convergence_table

# prepare the HTML template
tmplvars = {}
tmplvars['DATE'] = datetime.now().strftime("%c")
tmplvars['HOST'] = gethostname()
tmplvars['WHOAMI'] = getuser()

tmplvars['QUANTITY'] = quantity

# nice compact display of small and large floats, integers
compactfloat = lambda f: sub(r'\.0+', '',(u'%.3'+('f' if abs(f)<1e2 else 'e'))%f)

tmplvars['PARAMS_TABLE'] = paramtable.to_html()
tmplvars['ERROR_TABLE'] = errors.to_html()
tmplvars['CONVERGENCE_TABLE'] = convergence_table.to_html(float_format=compactfloat)

tmpl=open(report_templatefile, 'r').read().strip()
html = Template(tmpl).substitute(tmplvars)

with open(report_outputfile, 'w') as out:
	out.write(html)

print "Wrote report to %s." % report_outputfile

