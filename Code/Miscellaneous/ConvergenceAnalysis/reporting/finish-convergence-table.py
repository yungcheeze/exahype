#!/usr/bin/env python
#
# Try this script with
#   ipython -i finish-convergence-table.py
# to inspect all the lists and dicts.
#
# -- SK, 2016

import numpy as np
import pandas as pd
pd.set_option('expand_frame_repr', False) # display for debugging in terminal
# batteries:
import sys
from glob import glob
from os import path, stat
from re import match, sub
from os import path
import itertools, operator

# a module in this directory
from convergence_helpers import read_simulation_params, is_empty_file, simfile, gensvg, \
	Template, shortenPathsInTable, stripConstantColumns, keepColumnsIf, time_parser, RemoveStringInColumns

try:
	p = int(sys.argv[1])
	simulations = glob('simulations/p%d*/'%p) # mind the trailing slash to glob only directories
	report_outputfile = "simulations/generated-report-p%d.html"%p
except:
	print "Using all p orders"
	# determine finished simulations with
	# ./showSimulationProgress.sh  | grep finished | awk '{print $1}' > finished-simulations.txt
	simulations = [l.rstrip() for l in open('finished-simulations.txt')]
	report_outputfile = "simulations/generated-report.html"

# simulations is a list to dictionaries holding simulation data.

# which quantity shall we look at, in the moment?
quantity = 'error-bx.asc' # was error-rho.asc
quantityfile="output/" + quantity

# these are names for columns in our dataframes which can appear as columns
# in input and output data
idxSimName = 'SimulationName' # also used in showSimulationProgress.sh
idxParamFileName = 'ParamFileName'
idxQuantityFileName = 'QuantityFileName'

# What criterion do make to filter out bad simulations? These are examples:
minimalReductionsLengths = lambda SimRow: SimRow['EachRedLength'] > 120
emptyDirectoryCheck = lambda SimRow: not is_empty_file(SimRow['SimulationName'])
# now choose:
filterBadSimulations = minimalReductionsLengths

overview = Template("reporting/template-overview.html", report_outputfile).addStatistics()
evolution = Template("reporting/template-evolution.html", "simulations/evolution.html").addStatistics()
tpl = { 'QUANTITY': quantity }

print "Will work with %d simulations: " % len(simulations)
for t in enumerate(simulations): print " %i. %s" % t

## STEP 1: Retrieve Simulation parameters and statistics
## =====================================================

paramfiles = map(simfile('parameters.env'), simulations)
quantityfiles = map(simfile(quantityfile), simulations)

# 1A) In the first step, simtable will hold information about simulation paths
simtable = pd.DataFrame()
simtable[idxSimName] = simulations
simtable[idxParamFileName] = paramfiles
simtable[idxQuantityFileName] = quantityfiles
shortenPathsInTable(simtable, [idxParamFileName, idxQuantityFileName]) # for later html output

# 1B) read parameter files of each simulation
paramdicts = map(read_simulation_params, paramfiles) # list of dicts
paramtable = pd.DataFrame(paramdicts) # single table
paramtable[idxSimName] = simulations # give an index column

# 1C) Compute number of cells, etc.
porders = paramtable['EXAPORDER']
meshsizes = paramtable['EXAREALMESHSIZE'] # if not available, use 'EXAMESHSIZE'
widths = paramtable['EXASPEC_WIDTH']
# shall be all equal, of course.
if not np.all(widths == widths[0]):
	print "WARNING: Not all simulation domains are equal! Sizes are: ", widths
ncells = widths / meshsizes

# 1D) Beautifying of paramtable
# the following beautifying is only done for printing the paramtable
# beautify the parameter table: keep only columns with "EXA" inside
# paramtable = keepColumnsIf(paramtable, lambda c: 'EXA' in c)
# filter out PWD and OLDPWD which sometimes occur -.-
paramtable = keepColumnsIf(paramtable, lambda c: 'PWD' not in c)

# beautify the parameter table: shorten the paths
shortenPathsInTable(paramtable, ['EXABINARY', 'EXASPECFILE'])
# store number of cells back into paramtable for dumping
idxNcells = 'nCells'
paramtable[idxNcells] = ncells
tpl['SIMULATION_PARAMETERS_TABLE'] = paramtable.to_html()
# todo: Extract the "meaning" columns which hold documentation about the exa-columns.

# 1E) Read file about status of simulations (:= simulation statistics)
statisticstable = pd.read_csv("./logSimProgress.txt", sep="\t")
assert idxSimName in statisticstable.columns
# strip whitespace so a comparison is possible
statisticstable[idxSimName] = statisticstable[idxSimName].map(str.strip)
tpl['SIMULATION_STATISTICS_TABLE'] = statisticstable.to_html()
# todo: Create this on runtime without a shellscript

# 1F) Compose a shorter simulation table which contains paramtable and statistics
#     as well as shorten by extracting constant parameters
fullsimtable = pd.merge(simtable, pd.merge(paramtable, statisticstable, on=idxSimName), on=idxSimName)
reducedsimtable, constant_parameters = stripConstantColumns(fullsimtable)
# delete stuff we don't need, beautify up table
reducedsimtable = keepColumnsIf(reducedsimtable, lambda c: 'EXABINARY' not in c)
RemoveStringInColumns(reducedsimtable, 'EXA', inplace=True)

tpl['SIMULATION_TABLE'] = reducedsimtable.to_html()
tpl['CONSTANT_PARAMETERS'] = constant_parameters.to_html()

# 1G) Compute the total walltime
idxWalltime = 'Walltime' # as set in showSimulationProgress.sh
# standard time format (from /bin/time) is like "601m17.780s".
# We ignore the seconds and sum up the minutes
timedeltas = map(time_parser(r'\s*(?P<minutes>\d+)m'), fullsimtable[idxWalltime])
totaltime = reduce(operator.add, timedeltas)
totalhours = totaltime.total_seconds() / (60*60)
tpl['TOTAL_CPU_HOURS'] = ("%.1f" % totalhours)


## STEP 2: Load Error tables and compute convergence rate
## ======================================================

# 2A) Load the error table CSV files

# Caveats with header detection is very sensible to the first line's format.
# this will not work:
#    # 1:plotindex ,2:time ,3:l1norm ,4:l2norm ,5:max ,6:min ,7:avg ,
# in constrast, the line has to look like
#    1:plotindex 2:time 3:l1norm 4:l2norm 5:max 6:min 7:avg
# or even better
#    plotindex time l1norm l2norm max min avg
# which allows us to directly address the columns with their names.

# Read in of the actual error tables
errortables = [pd.read_csv(qf, delim_whitespace=True) for qf in quantityfiles]

# we do have one set of parameters for each error table
assert len(errortables) == len(paramtable)

# assign the index to each error table
idxcells = 'nCells' # name of 'number of cells' column
idxporder = 'pOrder' # name of 'polynomial order' column
for simncells,simporder,errortable in zip(ncells,porders,errortables):
	errortable[idxcells] = simncells
	errortable[idxporder] = simporder

# merge data frames to multiindex
errors = pd.concat(errortables).reset_index() #, keys=ncells, names=[idxcells])

# how much data rows can we compare? This is the number of maximum
# common evolution steps
simulationsPerNCells = errors.groupby(by=idxcells).size()
maxcomparisons = simulationsPerNCells.min()

print "For the following number of cells, we have this number of simlations:"
print simulationsPerNCells
print "So we can compare up to %d simulations for convergence analysis" % maxcomparisons

print "This is the full error table from all simulations (%d entries):" % len(errors)
print errors

idxplotindex = 'plotindex' # the column counting the rows in each simulation


# 2B) Compute the convergence rate

# we can either choose all overlapping data points for the common
# error (ceslice) or only the very last entry (celast). With the first
# one we can do plots showing the convergence order during evolution,
# with the last one we can show compact tables.
ceslice = errors.query('%s <= %d' % (idxplotindex, maxcomparisons))
celast = ceslice.groupby(by=[idxporder,idxcells], as_index=False).last().sort([idxporder,idxcells])
# either do groupby(..., as_index=False) and then .sort([idxporder,idxcells])
# or do as above: groupby() with index but still taking .last()
# celast: could also replace .last() by .tail(1) but would loose row index information
# CAVEAT: Make sure indices are unique! Otherwise the for ce.iterrows() will fail.
ce = ceslice

# @TODO: From here, go on with ce = celast for making
overview.set('COMBINED_FINAL_CONVERGENCE_ERRROR_TABLE', 'to be done with celast')
# @TODO: In the same time, go on with ceslice for making


# no index any more
# to get rid of the index:
#ce = ce.reset_index()
#ce.sort(idxcells, inplace=True) # inserts a "level_1" column

# find the reference nCells against which we make the convergence test.
# This is usually the one nCell at the same plotindex smaller than the
# actual one.
def findSmallerNcellsThan(num):
	smallerNcells = ncells[ ncells < num ]
	return np.nan if not len(smallerNcells) else smallerNcells.max()

idxprev = 'nSmaller' # column name of linked idxcells ('nCells')
ce[idxprev] = ce[idxcells].apply(findSmallerNcellsThan)

# insert the convergence measures for these columns
# we dropped 'min avg' as they don't tell us so much concerning convergence
columns = 'l1norm l2norm max'.split()
outcols = ["o"+c for c in columns]
idxtime = 'time' # the column for time

# add the new columns 
for newcol in outcols:
	ce[newcol] = 0.0

# Turn on or off row entanglement debugging
giveComparisonColumn = False
idxcomparisoncolumn="smallerRow"
if giveComparisonColumn:
	ce[idxcomparisoncolumn] = 'default'

# compute the actual convergence number for each column
# This requires unique indices in each row.
for rowindex, row in ce.iterrows():
	if (not row[idxprev]) or np.isnan(row[idxprev]):
		# there is no smaller resolution available
		if giveComparisonColumn:
			ce.set_value(rowindex, idxcomparisoncolumn, 'noSmaller')
		continue

	# this is the crucial search for the comparable row. We don't rely on
	# sorted rows but do an expensive search based on values here.
	targetrow = ce[(ce[idxcells] == row[idxprev]) & (ce[idxplotindex] == row[idxplotindex]) & (ce[idxporder] == row[idxporder])]
	if targetrow.empty:
		# could not find a previous step
		print "Did not find a target row for rowindex=%d and row[%s]=%s"  % (rowindex, idxprev, str(row[idxprev]))
		if giveComparisonColumn:
			ce.set_value(rowindex, idxcomparisoncolumn, 'CouldNotFind')
		continue

	if len(targetrow) != 1:
		print "Could not find unique counterpart for row ", row
		print "Found instead %d rows: " % len(targetrow), targetrow
		if giveComparisonColumn:
			ce.set_value(rowindex, idxcomparisoncolumn, 'FoundMultiple')
		continue

	# very verbose output
	#print "I am in row %d and matched nCells(%f) with prev(%f), outcome:" % (rowindex, row[idxcells], row[idxprev])
	#print pd.concat([ row.to_frame().T, targetrow ])

	for col, outcol in zip(columns, outcols):
		value = np.log(row[col] / targetrow[col]) / np.log( targetrow[idxcells] / row[idxcells] )
		#print "Computed row[%s]=%f" % (outcol, value)
		ce.set_value(rowindex, outcol, value)

	if giveComparisonColumn:
		ce.set_value(rowindex, idxcomparisoncolumn, str(list(targetrow[[idxplotindex,idxcells,idxtime,'l1norm']].values.flatten()))) # yes that's stupid

# print out a subset of the table
print "Comptuted this convergence table for the individual reductions"
print "(as l1norm, infnorm=max, etc.)"
convergence_table = ce[[idxporder,idxcells,idxprev,idxplotindex,idxtime] + columns + outcols + ([idxcomparisoncolumn] if giveComparisonColumn else [])]
convergence_table = convergence_table.sort([idxporder, idxcells, idxplotindex])
# do this for trying grouped output formatting:
#convergence_table = convergence_table.groupby(by=[idxporder, idxcells]).last()
#print convergence_table



# nice compact display of small and large floats, integers
#compactfloat = lambda f: sub(r'\.0+$', '',(u'%.3'+('f' if abs(f)<999 else 'e'))%f)
#compactfloat = lambda f: (u'%.3'+('f' if abs(f)<999 else 'e'))%f

# full tables which do *not* go to the overview

tpl['ERROR_EVOLUTION_TABLE'] = errors.to_html()
tpl['CONVERGENCE_EVOLUTION_TABLE'] = convergence_table.to_html() #float_format=compactfloat)





# we can generate plots, create strings holding the SVG file and embed
# the figures as inline SVG to the HTML file.

do_plots=True
if do_plots:
	print "Doing plots"
	import matplotlib
	matplotlib.use('Agg') # if headless
	import matplotlib.pyplot as plt
	plt.ion(); plt.clf()

	ycolumn = 'ol2norm'
	reprID = paramtable.at[0, 'EXAHYPE_INITIALDATA']
	reprSpecFile = paramtable.at[0, 'EXASPECFILE']

	fig = plt.figure(figsize=(18,8)) # plt.gcf()
	fig.suptitle("ExaHyPE convergence of %s for %s with %s" % (ycolumn, reprID, reprSpecFile), fontsize=18)

	uniquelist = lambda k: list(np.unique(k))
	lookupdict = lambda k,v: dict(zip(k,v))
	window = lambda l: itertools.izip(l, l[1:])
	iround = lambda f: int(round(f))

	allporders, allncells = uniquelist(porders), uniquelist(ncells)
	print "Plots are for porders=", allporders, " and ncells=", allncells

	possiblecolors = plt.cm.rainbow(np.linspace(0,1,len(allporders)))
	possiblestyles = ['-', '--', '-.', ':']
	possiblemarkers = ['o', 'v', 'x']
	# commented out because currently unused, but code works.
	#assert len(possiblemarkers) == len(uniquelist(ncells))
	#styles = lookupdict(allncells, possiblestyles)
	#markers = lookupdict(allncells, possiblemarkers)
	colors = lookupdict(allporders, possiblecolors)

	for Ni, (Nprev, Ncur) in enumerate(window(allncells)):
		plt.subplot(1, len(allncells)-1, Ni+1)

		plt.title("Convergence from %d to %d cells" % (iround(Nprev), iround(Ncur)))
		for (p,nc),rows in (convergence_table[convergence_table[idxcells]==Ncur]).sort([idxporder,idxcells], ascending=False).groupby(by=[idxporder,idxcells]):
			treshholdToShowPoints = 40
			style = "o-" if len(rows[ycolumn]) < treshholdToShowPoints else "-"				
			plt.plot(rows[idxtime], rows[ycolumn], style, label="P=%d" % int(p), color=colors[p])#, marker=markers[nc])

		plt.xlabel("Simulation time")
		plt.ylabel("Convergence order")
		plt.legend().draggable()
		plt.ylim(0,10)

	tpl['CONVERGENCE_SVG_FIGURE'] = gensvg(fig)
	tpl['ERROR_SVG_FIGURE'] = 'To be done'
else:
	tpl['CONVERGENCE_SVG_FIGURE'] = "<em>skipped plot generation</em>"
	tpl['ERROR_SVG_FIGURE'] = "<em>skipped plot generation</em>"

overview.execute(tpl, verbose=True)
evolution.execute(tpl, verbose=True)
