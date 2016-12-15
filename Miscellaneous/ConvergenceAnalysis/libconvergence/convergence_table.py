#!/usr/bin/env python
# -- SK, 2016

"""
This module is the key part in the libconvergence Python 2/3 library to
run exactness and convergence studies with ExaHyPE applications. It is
suitable to be run from command line or interactively from iPython.
"""

import numpy as np
import pandas as pd
pd.set_option('expand_frame_repr', False) # display for debugging in terminal
# batteries:
from glob import glob
from os import path, stat, getenv
from re import match, sub
import itertools, operator, sys, logging, inspect

# a module in this directory
from convergence_helpers import read_simulation_params, is_empty_file, simfile, gensvg, \
	Template, shortenPathsInTable, stripConstantColumns, keepColumnsIf, time_parser, \
	RemoveStringInColumns, is_headless, shell, StringIO
from convergence_arghelpers import ExaFrontend

def pdsort(df, by, ascending=True):
    "a pandas sort workaround which also works with pandas < 0.17"
    try:
        # pandas >= 0.17
        return df.sort_values(by=by, ascending=ascending)
    except:
        # pandas < 0.17
        return df.sort(by, ascending=ascending)

def merge_dicts(*dict_args):
    '''
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    '''
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result

def exitConvergenceStatus(convergencePassed):
	# other programs can use exit value to determine outcome of test
	sys.exit(0 if convergencePassed else -3)


class ConvergenceReporter:
	"""
	Does the actual convergence computation. Reads simulations from a list
	of succeeded simulations and uses the environment variable SIMBASE to
	determine the path to the simulations.
	"""

	def __init__(self):
		self.logger = logging.getLogger("ConvergenceReporter")

		# be_headless: do not show up matplotlib windows even if an X11
		# terminal is attached. Set False if you want to do work with plots.
		self.be_headless = True

		# run simulation statistics program before with
		# ./showSimulationProgress.sh  | grep FINISHED > simulations.txt
		# or so. Make sure it is a CSV table with column names, not only simulation name!
		self.simulationListFilename = False

		# add a path before the simulations names, if neccessary.
		self.simulationPathPrefix = getenv('SIMBASE', '')

		# which quantity shall we look at, in the moment?
		self.quantity = 'error-bx.asc' # was error-rho.asc

	def readSimulationProgress(self):
		"""
		Calls the shell script to read off the progress
		"""
		self.logger.info("Calling showSimulationProgress.sh to see progress")
		libconvergence_dir = path.dirname(__file__)
		script = path.join(libconvergence_dir, "showSimulationProgress.sh")
		csv = shell(script)
		self.logger.info("Simulation status:\n" + csv)
		return csv

	def add_group(self, argparser):
		group = argparser.add_argument_group('ConvergenceReporter', description=inspect.cleandoc(self.__doc__))
		group.add_argument('-o', '--porder', type=int, help="Do reporting only for one polynomial order. Default is all.")
		group.add_argument('-Q', '--quantity', type=str, default=self.quantity, help="Which quantity to look at")
		group.add_argument('-s', '--simulations', type=str, default=self.simulationListFilename, help="Path to CSV file "+
			"holding the simulations to consider. If not given, generates list by ./showSimulationProgress.sh on the fly")
		return group

	def apply_args(self, args, argparser):
		if args.porder:
			p = args.porder
			self.logger.info("Creating convergence table only for p=%d" % p)
			self.simulations = glob('simulations/p%d*/'%p) # mind the trailing slash to glob only directories
			self.report_outputfile = "simulations/generated-report-p%d.html"%p
		else:
			self.logger.info("Using all p orders")
			self.simulations = None
			self.report_outputfile = "simulations/generated-report.html"

		if args.simulations:
			self.logger.info("Reading off simulations from file %s" % args.simulations)
			self.simulationListFilename = args.simulations
		else:
			self.simulationListFilename = StringIO(self.readSimulationProgress())

		self.quantity = args.quantity

	def start(self):
		# copy variables in local namespace for ... yeah.
		simulations = self.simulations

		# simulations is a list to dictionaries holding simulation data.
		self.quantityfile = "output/" + self.quantity

		# these are names for columns in our dataframes which can appear as columns
		# in input and output data
		idxSimName = 'SimulationName' # also used in showSimulationProgress.sh
		idxParamFileName = 'ParamFileName'
		idxQuantityFileName = 'QuantityFileName'

		# What criterion do make to filter out bad simulations? These are examples:
		minimalReductionsLengths = lambda SimRow: SimRow['EachRedLength'] > 120
		emptyDirectoryCheck = lambda SimRow: not is_empty_file(SimRow['SimulationName'])
		# now choose:
		goodSimulations = minimalReductionsLengths

		# Load table of simulations
		# The shellscript takes around ~ 10-60 Seconds to generate the data, therefore
		# the generation is offloaded.
		statisticstable = pd.read_csv(self.simulationListFilename, sep="\t")
		assert idxSimName in statisticstable.columns, "statisticstable misses essential data: "+str(statisticstable)

		# strip whitespace so a comparison is possible
		statisticstable[idxSimName] = statisticstable[idxSimName].map(str.strip)

		# allow to join a path before if simnames are somewhat broken or so.
		if self.simulationPathPrefix:
		    self.logger.info("Applying prefix '%s' on each simulation name" % self.simulationPathPrefix)
		    prefixer = lambda simname: path.join(self.simulationPathPrefix, simname)
		    # if working with pandas:
		    statisticstable[idxSimName] = statisticstable[idxSimName].apply(prefixer)
		    # if working with the list, interchange with if not simulations: ... but doesn't work.
		    #simulations = map(prefixer, simulations)

		if not simulations:
			simulations = list(statisticstable[idxSimName])

		overview = Template("reporting/template-overview.html", self.report_outputfile).addStatistics()
		evolution = Template("reporting/template-evolution.html", "simulations/evolution.html").addStatistics()
		overview.set('LINK_DETAILED_REPORT', path.basename(evolution.outputfile)) # link them together
		tpl = { 'QUANTITY': self.quantity } # common template variables

		self.logger.info("Will work with %d simulations: " % len(simulations))
		for t in enumerate(simulations): self.logger.info(" %i. %s" % t)

		## STEP 1: Retrieve Simulation parameters and statistics
		## =====================================================

		paramfiles = map(simfile('parameters.env'), simulations)
		quantityfiles = map(simfile(quantityfile), simulations)

		# 1A) In the first step, simtable will hold information about simulation paths
		simtable = pd.DataFrame()
		simtable[idxSimName] = simulations
		simtable[idxParamFileName] = paramfiles
		simtable[idxQuantityFileName] = quantityfiles

		# 1B) read parameter files of each simulation
		paramdicts = map(read_simulation_params, paramfiles) # list of dicts
		paramtable = pd.DataFrame(paramdicts) # single table
		paramtable[idxSimName] = simulations # give an index column

		# 1C) Compute number of cells, etc.

		def case_insensitive(df, colname):
		    "find a column name by case-insensitive matching"
		    col_list = list(df)
		    try:
			# this uses a generator to find the index if it matches, will raise an exception if not found
			return col_list[next(i for i,v in enumerate(col_list) if v.lower() == colname.lower())]
		    except:
			raise KeyError("Could not find column '%s' in list of columns '%s'" % (colname, str(col_list)) )

		# A case-insensitive finder for paramtable
		ci_paramtable = lambda colname: case_insensitive(paramtable, colname)
		get_ci_paramtable = lambda colname: paramtable[ci_paramtable(colname)]

		porders = get_ci_paramtable('EXAPORDER')
		meshsizes = get_ci_paramtable('EXAREALMESHSIZE') # if not available, use 'EXAMESHSIZE'
		widths = get_ci_paramtable('ExaWidth') # used to be calles EXASPEC_WIDTH
		# shall be all equal, of course.
		if not np.all(widths == widths[0]):
			self.logger.warning("Not all simulation domains are equal! Sizes are: ", widths)
		ncells = widths / meshsizes

		# 1D) Beautifying of paramtable
		# the following beautifying is only done for printing the paramtable
		# beautify the parameter table: keep only columns with "EXA" inside
		## paramtable = keepColumnsIf(paramtable, lambda c: 'EXA' in c)
		# filter out PWD and OLDPWD which sometimes occur -.-
		paramtable = keepColumnsIf(paramtable, lambda c: 'PWD' not in c)

		# beautify the parameter table: shorten the paths
		shortenPathsInTable(paramtable, map(ci_paramtable, ['EXABINARY', 'EXASPECFILE']))
		linkifyFormatter = { idxSimName: lambda l: u"<a href='%s'>%s</a>"%(l,l) }
		defaultFormatter = { col: None for col in paramtable.columns.values }
		# store number of cells back into paramtable for dumping
		idxNcells = 'nCells'
		paramtable[idxNcells] = ncells
		tpl['SIMULATION_PARAMETERS_TABLE'] = paramtable.to_html(formatters=merge_dicts(defaultFormatter, linkifyFormatter))
		tpl['SIMULATION_STATISTICS_TABLE'] = statisticstable.to_html()

		# todo: Extract the "meaning" columns which hold documentation about the exa-columns.

		# 1F) Compose a shorter simulation table which contains paramtable and statistics
		#     as well as shorten by extracting constant parameters
		fullsimtable = pd.merge(simtable, pd.merge(paramtable, statisticstable, on=idxSimName), on=idxSimName)

		# Compute the filter of entries which will be processed in the further steps
		includedSimulations = fullsimtable.apply(goodSimulations, axis=1) # row-wise
		idxIgnoredSims='IgnoredSimulations'
		IgnoredKeyworder = lambda b: "Included" if b else "Discarded"
		fullsimtable[idxIgnoredSims] = includedSimulations.map(IgnoredKeyworder)

		# shorten paths for reduction
		shortenPathsInTable(fullsimtable, [idxParamFileName, idxQuantityFileName])
		reducedsimtable, constant_parameters = stripConstantColumns(fullsimtable)
		# delete stuff we don't need, beautify up table
		reducedsimtable = keepColumnsIf(reducedsimtable, lambda c: 'EXABINARY' not in c)
		RemoveStringInColumns(reducedsimtable, 'EXA', inplace=True)

		tpl['SIMULATION_TABLE'] = reducedsimtable.to_html()
		tpl['CONSTANT_PARAMETERS'] = constant_parameters.to_html()

		# 1H) Compute the total walltime
		idxWalltime = 'Walltime' # as set in showSimulationProgress.sh
		# standard time format (from /bin/time) is like "601m17.780s".
		# We ignore the seconds and sum up the minutes
		timedeltas = map(time_parser(r'\s*(?P<minutes>\d+)m'), fullsimtable[idxWalltime])
		badlyparsedCriterion = lambda timedelta: not timedelta
		correctlyparsed = map(badlyparsedCriterion, timedeltas)
		errnousentries = len(timedeltas) - len(correctlyparsed) # find not correctly parsed entries
		timedeltas = correctlyparsed
		totaltime = reduce(operator.add, timedeltas)
		try:
		    totalhours = totaltime.total_seconds() / (60*60)
		    tpl['TOTAL_CPU_HOURS'] = ("%.1f" % totalhours) + (" (ignoring %i errnous entries)"%errnousentries if errnousentries else "")
		except:
		    tpl['TOTAL_CPU_HOURS'] = 'Not determinable'


		## STEP 2: Load Error tables and compute convergence rate
		## ======================================================

		# 2A) Filter out all bad entries which are no further interesting for automatic
		#     processing
		paramtable = paramtable[includedSimulations]
		quantityfiles = simtable[includedSimulations][idxQuantityFileName]

		# 2B) Load the error table CSV files

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

		self.logger.info("For the following number of cells, we have this number of simlations:")
		self.logger.info(simulationsPerNCells)
		self.logger.info("So we can compare up to %d simulations for convergence analysis" % maxcomparisons)

		# too much data:
		#print "This is the full error table from all simulations (%d entries):" % len(errors)
		#print errors

		idxplotindex = 'plotindex' # the column counting the rows in each simulation


		# 2B) Compute the convergence rate

		# we can either choose all overlapping data points for the common
		# error (ceslice) or only the very last entry (celast). With the first
		# one we can do plots showing the convergence order during evolution,
		# with the last one we can show compact tables.
		ceslice = errors.query('%s <= %d' % (idxplotindex, maxcomparisons))
		celast = pdsort(ceslice.groupby(by=[idxporder,idxcells], as_index=False).last(), by=[idxporder,idxcells])
		# either do groupby(..., as_index=False) and then .sort([idxporder,idxcells])
		# or do as above: groupby() with index but still taking .last()
		# celast: could also replace .last() by .tail(1) but would loose row index information
		# CAVEAT: Make sure indices are unique! Otherwise the `for row in ce.iterrows()` will fail.
		ce = ceslice

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

		#def compute_convergence_order(ce):
		#	"ce: common error table"
		#	ce = ce.copy(deep=True)
		idxprev = 'nSmaller' # column name of linked idxcells ('nCells')
		ce[idxprev] = ce[idxcells].apply(findSmallerNcellsThan)

		# insert the convergence measures for these columns
		# we dropped 'min avg' as they don't tell us so much concerning convergence
		errorColumns = 'l1norm l2norm max'.split()
		rateColumns = ["o"+c for c in errorColumns]
		idxtime = 'time' # the column for time

		# add the new columns 
		for newcol in rateColumns:
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
				self.logger.error("Did not find a target row for rowindex=%d and row[%s]=%s"  % (rowindex, idxprev, str(row[idxprev])))
				if giveComparisonColumn:
					ce.set_value(rowindex, idxcomparisoncolumn, 'CouldNotFind')
				continue

			if len(targetrow) != 1:
				self.logger.error("Could not find unique counterpart for row ", row)
				self.logger.error("Found instead %d rows: " % len(targetrow), targetrow)
				if giveComparisonColumn:
					ce.set_value(rowindex, idxcomparisoncolumn, 'FoundMultiple')
				continue

			# very verbose output
			#print "I am in row %d and matched nCells(%f) with prev(%f), outcome:" % (rowindex, row[idxcells], row[idxprev])
			#print pd.concat([ row.to_frame().T, targetrow ])

			for col, outcol in zip(errorColumns, rateColumns):
				value = np.log(row[col] / targetrow[col]) / np.log( targetrow[idxcells] / row[idxcells] )
				#print "Computed row[%s]=%f" % (outcol, value)
				ce.set_value(rowindex, outcol, value)

			if giveComparisonColumn:
				ce.set_value(rowindex, idxcomparisoncolumn, str(list(targetrow[[idxplotindex,idxcells,idxtime,'l1norm']].values.flatten()))) # yes that's stupid


		self.logger.info("Computing convergence tables...")
		#convEvolution = compute_convergence_order(ceslice)
		#convFinal = compute_convergence_order(celast)

		# nice compact display of small and large floats, integers
		#compactfloat = lambda f: sub(r'\.0+$', '',(u'%.3'+('f' if abs(f)<999 else 'e'))%f)
		#compactfloat = lambda f: (u'%.3'+('f' if abs(f)<999 else 'e'))%f

		# full tables which do *not* go to the overview

		# print out a subset of the table
		self.logger.info("Comptuted this convergence table for the individual reductions")
		self.logger.info("(as l1norm, infnorm=max, etc.)")
		convergence_table = ce[[idxporder,idxcells,idxprev,idxplotindex,idxtime] + errorColumns + rateColumns + ([idxcomparisoncolumn] if giveComparisonColumn else [])]
		self.logger.info( str(convergence_table) )
		convergence_table = pdsort(convergence_table, by=[idxporder, idxcells, idxplotindex])
		final_time_convergence_table = convergence_table.groupby(by=[idxporder, idxcells]).last()

		# We can also programmatically decide whether decent convergence
		# is present or not. We measure acceptability as the difference from
		# ideal scaling.
		def acceptability(row):
			# assuming the ideal scaling is always bigger than the real scaling, we sum
			# up defects which are an indicator about scaling quality
			return sum([ (row[idxporder]+1) - row[rateCol] for rateCol in rateColumns ])

		convergenceQuality = convergence_table.apply(acceptability, axis=1)
		convergenceQualityIndicator = convergenceQuality.mean()

		# todo: We can give names
		convergenceQualityBelowNames = { 1: 'excellent', 5: 'passes', 10: 'failed' }
		convergencePassed = ( convergenceQualityIndicator < 5.0 )

		tpl['CONVERGENCE_QUALITY_INDICATOR'] = "%.2f" % convergenceQualityIndicator
		tpl['CONVERGENCE_PASSED'] = ('PASSED' if convergencePassed else 'FAILED')

		tpl['ERROR_EVOLUTION_TABLE'] = errors.to_html()
		tpl['CONVERGENCE_EVOLUTION_TABLE'] = convergence_table.to_html() #float_format=compactfloat)
		tpl['COMBINED_FINAL_CONVERGENCE_ERRROR_TABLE'] = final_time_convergence_table.to_html()


		# we can generate plots, create strings holding the SVG file and embed
		# the figures as inline SVG to the HTML file.

		do_plots=True
		if do_plots:
			self.logger.info("Doing plots")
			tpl.update(self.do_plots())
		else:
			self.logger.info("Skipped plot creation")
			tpl['CONVERGENCE_SVG_FIGURE'] = "<em>skipped plot generation</em>"
			tpl['ERROR_SVG_FIGURE'] = "<em>skipped plot generation</em>"

		overview.execute(tpl, verbose=True)
		evolution.execute(tpl, verbose=True)

		self.logger.info("Convergence test results:")
		self.logger.info("Convergence factor is %.2f" % convergenceQualityIndicator)
		self.logger.info("Convergence test is " + ("PASSED" if convergencePassed else "FAILED"))

		return convergencePassed

	def do_plots(self):
		import matplotlib
		if self.be_headless or is_headless():
			matplotlib.use('Agg')
		import matplotlib.pyplot as plt
		plt.ion(); plt.clf()

		uniquelist = lambda k: list(np.unique(k))
		uniqueintlist = lambda k: map(int, uniquelist(k))
		lookupdict = lambda k,v: dict(zip(k,v))
		window = lambda l: itertools.izip(l, l[1:])
		iround = lambda f: int(round(f))
		# repeat tile until it fits to target. or reduce tile down to target.
		repeatfit = lambda tile, target: (max(1,int(len(target)/len(tile)))*tile)[:len(target)]
		# apply repeatfit to the lookup values of lookupdict
		lookupfitdict = lambda k,v: lookupdict(k, repeatfit(v, k))

		allporders, allncells = uniqueintlist(porders), uniqueintlist(ncells)
		self.logger.info("Plots are for porders=%s and ncells=%s" % (allporders,allncells))

		# this never gave good colors:
		#possiblecolors = plt.cm.jet(np.linspace(0,0.9,len(allporders)))
		# use instead:
		possiblecolors = "r g k b m y r c m".split(" ")
		#possiblestyles = ['-', '--', '-.', ':']
		possiblemarkers = ['o', '>', 's', 'D','H', '*','<',]
		# commented out because currently unused, but code works.
		#assert len(possiblemarkers) == len(uniquelist(ncells))
		#styles = lookupdict(allncells, possiblestyles)
		markerdict = lookupfitdict(allporders, possiblemarkers)
		colordict = lookupfitdict(allncells, possiblecolors)
		markers = lambda p, nc: markerdict[int(p)]
		colors = lambda p, nc: colordict[int(nc)]

		##### SIMPLE ERROR EVOLUTION PLOTS
		errorPlot = plt.figure(figsize=(18,8))
		ycolumn = 'max' # l2norm, l1norm

		num_markers = 10
		dt_markers = max(errors[idxtime]) / num_markers
		for groupidx,rows in errors.groupby(by=[idxporder,idxcells]):
			# groupidx = (p,nc)
			plt.plot(rows[idxtime], rows[ycolumn], label="P=%d,Nc=%d"%groupidx, color=colors(*groupidx), marker=markers(*groupidx), markevery=0.03)

		plt.legend(loc='center left', bbox_to_anchor=(1,0.5))
		plt.subplots_adjust(right=0.8)
		plt.title("Error evolution: %s of quantity %s" % (ycolumn, self.quantity))
		plt.xlabel("Simulation time")
		plt.ylabel("Error")
		ax = plt.gca()
		ax.set_yscale('log')
		plt.ylim(1e-10, 1e0)

		##### CONVERGENCE EVOLUTION PLOTS

		ycolumn = 'ol2norm'
		reprID = get_ci_paramtable('EXAHYPE_INITIALDATA').iloc[0]
		reprSpecFile = get_ci_paramtable('EXASPECFILE').iloc[0]

		convergencePlot = plt.figure(figsize=(18,8)) # plt.gcf()
		convergencePlot.suptitle("ExaHyPE convergence of %s for %s with %s" % (ycolumn, reprID, reprSpecFile), fontsize=18)

		# this is an ugly bug, we are dealing with ints but storing floats
		# and then making comparisons float(9) == float(9.000000...) which never hold true.
		allporders, allncells = uniquelist(porders), uniquelist(ncells)

		# do plots for all neighbouring (Nprev, Ncur)
		for Ni, (Nprev, Ncur) in enumerate(window(allncells)):
			plt.subplot(1, len(allncells)-1, Ni+1)

			plt.title("Convergence from %d to %d cells" % (iround(Nprev), iround(Ncur)))
			for groupidx,rows in pdsort(convergence_table[convergence_table[idxcells]==Ncur], by=[idxporder,idxcells], ascending=False).groupby(by=[idxporder,idxcells]):
				treshholdToShowPoints = 40
				style = "o-" if len(rows[ycolumn]) < treshholdToShowPoints else "-"				
				(p,nc) = groupidx
				plt.plot(rows[idxtime], rows[ycolumn], style, label="P=%d" % int(p))#, color=colors(*groupidx))#, marker=markers[nc], markersize=25)

			plt.xlabel("Simulation time")
			plt.ylabel("Convergence order")
			plt.legend()#.draggable()
			plt.ylim(0,10)


		tpl = {}
		try:
			tpl['CONVERGENCE_SVG_FIGURE'] = gensvg(convergencePlot)
			tpl['ERROR_SVG_FIGURE'] = gensvg(errorPlot)
		except Exception as e:
			self.logger.error("Could not create plots, due to exception: %s" % str(e))
			tpl['CONVERGENCE_SVG_FIGURE'] = "<em>Creation crashed!</em>"
			tpl['ERROR_SVG_FIGURE'] = "<em>Creation crashed!</em>"
		return tpl

def main():
	"A main function for the convergence table computation."
	logger = logging.getLogger("convergence_table.py")

	frontend = ExaFrontend(program_description=__doc__)
	reporter = ConvergenceReporter()
	frontend.add_module(reporter)

	args = frontend.parse_args()
	convergencePassed = reporter.start()

	logger.info("Finishing with convergencePassed=%s" % str(convergencePassed))
	exitConvergenceStatus(convergencePassed)

if __name__ == "__main__":
	main()


