#!/usr/bin/env python
# Run convergence tests.
#
# This is basically python version agnostic and should run with both python2
# and python3.
#
# ExaHYPE 2016, SvenK

"""
The Convergence starter package contains a number of classes to generally
manage ExaHyPE runs driven by environmental variables and templated spec files
and especially to manage polynomial convergence tests.
All those classes are abstract and supposed to be used with a concrete frontend
python program.
"""

from numpy import array
import os, sys, logging # batteries

# libconvergence, same directory:
from convergence_table import ConvergenceReporter, exitConvergenceStatus
from convergence_helpers import shell, getenv, pidlist, runBinary, cleandoc, Future

class ConvergenceTest:
	"""
	A convergence simulation class basically holds a "settings" dict
	which becomes the environment variables for an invocation of the
	runTemplatedSpecfile.sh shell script.

	We only need this shell script as it is quite convenient to setup
	the directory layout, etc. from a shell script.
	"""

	def __init__(self, name="UnnamedTest"):
		self.name = name
		self.logger = logging.getLogger(repr(self))
		self.boundsettings = {}

		settings = {}

		settings['SIMBASE'] = getenv('SIMBASE', default="simulations/")

		# a crude way to determine the base directory
		exabase = settings['ExaBase'] = shell("git rev-parse --show-toplevel")
		examisc = settings['ExaMisc'] = os.path.join(exabase, "Miscellaneous/")

		settings['ExaRunner'] = os.path.join(examisc, "RunScripts/runTemplatedSpecfile.sh")

		# template to set up a queueing system
		# an example is:
		#   "srun -n1 --partition=x-men --time=29:00:00 --mem=0 --job-name=p{ExapOrder}-m{ExaMeshSize}-ConvergenceStudies-Euler"
		settings['QRUNTPL'] = getenv('QRUNTPL', default='')

		self.settings = settings

	def stringify(self, inplace=False):
		"Ensure self.settings is suitable for being an ENV array"
		stringified = { k: str(v) for k,v in self.settings.iteritems() }
		if inplace:
			self.settings = stringified
		return stringified

	def __repr__(self):
		return "ConvergenceTest(%s)" % self.name

	def add_group(self, parser):
		"""
		Overwrite this function to add argument parsers or groups, like in
		https://docs.python.org/2/library/argparse.html#argument-groups
		See the PolyorderTest as an example.
		"""
		pass


	def apply_args(self, args, argparser):
		"""
		Overwrite this function for an early check of the parameters
		"""
		pass

	def start(self):
		"""
		Starts a simulation. This method should be basically overloaded
		with something useful replacing some parameters, etc.
		"""
		self.check()
		env = os.environ.copy()
		env.update(self.settings)
		
		self.logger.info("Generic ConvergenceTest starts {ExaRunner}".format(**env))
		return runBinary(env['ExaRunner'], env)
	
	def startFromArgs(self, args, parser):
		raise NotImplementedError("Overwrite this function to start the test from ConvergenceFrontend.")

class ParametricTest(ConvergenceTest):
	"""
	Allows to pass command line options to settings in the ConvergenceTest.
	
	Example usage:
	
	> test = ParametricTest("DemoTest")
	> test.settings['ExaStaticFoo'] = 3.141
	> test.settings['ExaFrequentlyChangingBar'] = test.commandline('-b', '--bar', choices={ 1,2,3 }, help="Typical Bar values")
	> test.settings['ExaAlsoChangingBaz'] = test.commandline('-z', '--baz', type=str, help="Just another argument")
	
	This mechanism breaks with direct starting from ConvergenceTest:
	> test.start()
	
	However, it goes nicely with the Frontend and command line parsing:
	
	> ConvergenceFrontend(test, description="My Demonstrator")
	
	Side notice: Somehow this is really overengineered.
	"""
	futureclass = Future
	
	def __init__(self, name="UnnamedParametricTest"):
		ConvergenceTest.__init__(self, name)
	
	def commandline(self, *args, **kwargs):
		"""
		Neat syntactic sugar to return an instance of Future which is later used
		to call the bound method of argparser.ArgumentParser.add_argument.
		"""
		return self.futureclass(*args, **kwargs)

	def add_group(self, parser):
		"Hook the future settings into the parser"
		group = parser.add_argument_group(title=repr(self), description=cleandoc(self))
		for key,future in self.settings.iteritems():
			if isinstance(future, self.futureclass):
				future.apply(group.add_argument)
		return group

	def apply_args(self, args, argparser):
		"Reads out the argument line values for the future settings and pass them"
		argns = vars(args) # namespace->dict
		for key,future in self.settings.iteritems():
			if isinstance(future, self.futureclass):
				# extract from the argparser namespace
				self.settings[key] = argns[future.value().dest]

class PolyorderTest(ParametricTest):
	"""
	The PolyorderTest is a convergence test on a uniform grid where both
	the grid spacing as well as the polynomial order in ADERDG is changed
	in order to see convergence accross both degrees of freedom.
	"""
	def __init__(self, name="UnnamedPolyorderTest", polyorders=[], width=0.0, depths=[]):
		ConvergenceTest.__init__(self, name)

		# typical values:
		#polyorders = arange(2,10)
		#width = 2.
		#depths = arange(1,6)

		self.polyorders = polyorders
		self.width = width
		self.depths = depths

		self.meshsize = lambda width, depth: width / 3.**depth  # actual meshsize we will get (dx)
		self.numcells = lambda width, meshsize: width / meshsize
		self.maxmeshsizefactor = 1.1 # to ensure meshsize is little bit smaller
		self.maxmeshsize = lambda meshsize: meshsize * self.maxmeshsizefactor

		# defaults paths
		self.settings['ExaBinaryTpl'] = '{ExaBinary}-p{ExapOrder}'
		self.settings['SIMDIRTpl'] = "{SIMBASE}/p{ExapOrder}-meshsize{ExaMeshSize}/"
		
		# future data binding.
		# "meshsize" and "polyorder" as argsNS members are also directly used in startFromArgs
		self.settings["ExaRealMeshSize"] = self.commandline('-m', '--meshsize', type=float, help="Maximum mesh size (dx)")
		self.settings['ExapOrder'] = self.commandline('-p', '--polyorder', type=int, help="Polynomial order")
		# Since I created this ParametricTest, now I have to deal with it...

	def until(self, maxcells):
		"Small utility function" 
		return self.res[ self.res['numcells'] <= maxcells ]

	def computeCombinations(self):
		import pandas as pd
		self.res = pd.DataFrame({'depth': array(self.depths)})
		self.res['meshsize'] = dx = self.meshsize(self.width, array(self.depths))
		self.res['numcells'] = self.numcells(self.width, dx)
		self.res['maxmeshsize'] = self.maxmeshsize(dx)

		# default runrange:
		until = self.until
		# adapt the number of runs/ maximum cells for the polynomial degree
		# adaptrunrange[<polyorder>] -> [<subset of maxmeshsizes>]
		# numcells is array([   3.,    9.,   27.,   81.,  243.])
		self.adaptiveRunRange = { 2: until(243), 3: until(243), 4: until(243), 5: until(243),
				6: until(81), 7: until(81),
				8: until(27), 9: until(27) }

		return self.res

	def start(self, polyorder, maxmeshsize, meshsize):
		env = self.stringify()
		env['ExaMeshSize'] = str(maxmeshsize)
		env["ExaRealMeshSize"] = str(meshsize)
		env['ExapOrder'] = str(polyorder)

		env['ExaBinary'] = env['ExaBinaryTpl'].format(**env)
		env['SIMDIR']= env['SIMDIRTpl'].format(**env)

		if 'QRUNTPL' in env:
			env['QRUN'] = env['QRUNTPL'].format(**env)
		
		#if not os.path.exists(env['ExaBinary']):
		#		raise IOError("Failure: '{ExaBinary}' does not exist, probably you forgot to"
		#			"compile for all neccessary polynomial orders before?"
		#			"Hint: Try `exa polycompile EulerFlow`".format(**env))


		env['BaseNameExaRunner'] = os.path.basename(env['ExaRunner'])
		self.logger.info("Starting p={ExapOrder}, maxmeshsize={ExaMeshSize} with {BaseNameExaRunner}".format(**env))
		return runBinary(env['ExaRunner'], env)


	def startRange(self, adaptrunrange=None):
		"""
		Having polyorders, depths being two lists, run the cartesian product
		itertools.product(polyorders, depths) of simulations.
		However, do not run too fine simulations: Make a cut off at high polynomial
		orders.
		If you call this with no arguments, it takes the current value of self.runrange.
		If you have not set self.runrange by yourself, it defaults to a reasonable set
		of simulations.
		"""
		self.logger.info("ExaHyPE Polynomial Convergence Analysis")
		self.logger.info("=======================================")

		self.logger.info("Can do convergence analysis on %fx%f sized domain with the following grid props: " % (self.width,self.width))
		self.logger.info("\n" + str(self.res))
		self.logger.info("Will do all these tests for these orders of the polynomial order:")
		self.logger.info(self.polyorders)
		self.logger.info("However, now reducing the runs at big resolutions")
	
		processes = []
		if not adaptrunrange:
			adaptrunrange = self.adaptiveRunRange
		for p, rows in adaptrunrange.iteritems():
			for i, row in rows.iterrows():
				# for testing:
				#proc = subprocess.Popen(["/bin/sleep", str(p)])
				# for real:
				proc = self.start(p, row['maxmeshsize'], row['meshsize'])
				processes.append(proc)
		return processes

	def add_group(self, parser):
		group = ParametricTest.add_group(self, parser)
		# now implemented as futures:
		#group.add_argument('-p', '--polyorder', type=int, help="Polynomial order")
		#group.add_argument('-m', '--meshsize', type=float, help="Maximum mesh size (dx)")
		group.add_argument('-a', '--all', action='store_true', help='Start all sensible simulations')

	def startFromArgs(self, args, parser):
		if args.polyorder or args.meshsize:
			if(not args.polyorder or not args.meshsize):
				parser.error("Please specify both the polynomial order and the meshsize or none")
			p = self.start(args.polyorder, self.maxmeshsize(args.meshsize), args.meshsize)
			self.logger.info("ExaHyPE runner has been started, PID = %d" % p.pid)
			return [p]
		elif args.all:
			processes = self.startRange()
			self.logger.info("%d processes have been started with PIDS: " % len(processes))
			self.logger.info(pidlist(processes))
			return processes
		#elif args.nostart:
		#	self.logger.info("Skipping all ExaHyPE process starting.")
		#	return []
		else:
			self.logger.error("Choose either --all, or -p/-m combination for run")
			parser.print_help()
			sys.exit(2)



