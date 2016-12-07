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

from numpy import arange, array
import subprocess, os, sys, time, logging, inspect # batteries
logger = logging.getLogger("convergence_starter")

# same directory helpers
from convergence_arghelpers import ExaFrontend
#import convergence_table

# helpers:
shell = lambda cmd: subprocess.check_output(cmd, shell=True).strip()
getenv = lambda key, default=None: os.environ[key] if key in os.environ else default
#me = os.path.basename(__file__)
pidlist = lambda processes: " ".join([str(proc.pid) for proc in processes])

def runBinary(binary, envUpdate):
	"""
	Runs the command "binary", this may be only one program without parameters,
	with the environment which is composed of the environmental variables and
	the envUpdate dictionary.
	Example usage:
	> runBinary("ExaHyPE-FooBar", { 'EXAHYPE_SKIP_TEST': 'True' })

	You can replace subprocess.Popen with other alternatives, but make sure
	they return some kind of process handle which is subsequently used
	"""
	env = os.environ.copy()
	env.update(envUpdate)
	
	if not os.path.exists(binary):
		raise IOError("Failure: '%s' does not exist" % binary)

	return subprocess.Popen([binary], env=env)

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

class PolyorderTest(ConvergenceTest):
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
		self.runrange = { 2: until(243), 3: until(243), 4: until(243), 5: until(243),
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
		logger.info("ExaHyPE Polynomial Convergence Analysis")
		logger.info("=======================================")

		logger.info("Can do convergence analysis on %fx%f sized domain with the following grid props: " % (self.width,self.width))
		logger.info(self.res)
		logger.info("Will do all these tests for these orders of the polynomial order:")
		logger.info(self.polyorders)
		logger.info("However, now reducing the runs at big resolutions")
	
		processes = []
		if not adaptrunrange:
			adaptrunrange = self.adaptrunrange
		for p, rows in adaptrunrange.iteritems():
			for i, row in rows.iterrows():
				# for testing:
				#proc = subprocess.Popen(["/bin/sleep", str(p)])
				# for real:
				proc = self.start(p, row['maxmeshsize'], row['meshsize'])
				processes.append(proc)
		return processes

	def add_group(self, parser):
		group = parser.add_argument_group(title=repr(self), description=inspect.cleandoc(self.__doc__))
		group.add_argument('-p', '--polyorder', type=int, help="Polynomial order")
		group.add_argument('-m', '--meshsize', type=float, help="Maximum mesh size (dx)")
		group.add_argument('-a', '--all', action='store_true', help='Start all sensible simulations')

	def startFromArgs(self, args, parser):
		if args.polyorder or args.meshsize:
			if(not args.polyorder or not args.meshsize):
				parser.error("Please specify both the polynomial order and the meshsize or none")
			p = self.start(args.polyorder, self.maxmeshsize(args.meshsize), args.meshsize)
			self.logger.info("ExaHyPE runner has been started, PID = %d" % p.pid)
			return [p]
		elif args.all:
			processes = runRange()
			self.logger.info("%d processes have been started with PIDS: " % len(processes))
			self.logger.info(pidlist(processes))
			return processes
		else:
			self.logger.error("Choose either --all or -p/-m combination for run")
			parser.print_help()
			sys.exit(2)

class ConvergenceStarterReportingAdapter:
	"""
	This class allows direct invocation of the reporting by waiting
	for the convergence studies to be finished.
	The reporting can always also be called manually by invoking the
	CLI program.
	"""

	def __init__(self):
		self.logger = logging.getLogger("ConvergenceStarterReportingAdapter")

	def add_group(self, parser):
		group = parser.add_argument_group(title="Reporting", description=inspect.cleandoc(self.__doc__))
		group.add_argument('-w', '--wait', action='store_true', help="Wait until ExaHyPE processes have finished")
		group.add_argument('-r', '--reporting', action='store_true', help="Do the reporting. Will trigger -w.")

	def apply_args(self, args, argparser):
		self.reporting = args.reporting
		self.wait = args.wait

		if args.reporting and not args.wait:
			self.logger.info("Switching on --wait")
			self.wait = True

	def dispatchProcesses(self, processes):
		if self.wait:
			livingprocs = processes[:] # copy
			while True:
				# filter out finished processes
				livingprocs = [proc for proc in livingprocs if proc.poll() == None]
				self.logger.info("Waiting for %d processes (%s) to finish..." % (len(livingprocs),pidlist(livingprocs)))
				if not len(livingprocs):
					break
				time.sleep(1)
			
				# Todo: Should kill processes when this script is killed inside this loop.
				# or tell processes to quit when args.wait is on and parent is killed.
			
				# do this to access values
				#proc.communicate()

			exitcodes = [proc.returncode for proc in processes]
			exitcode = sum(exitcodes)
			self.logger.info("All processes finished with summed %d. Exit codes are: %s" % (exitcode, str(exitcodes)))

			if self.reporting:
				self.convergencePassed = self.startConvergenceTable()
		else:
			self.logger.info("Convergence tests running in background, immediately finishing.")

	def startConvergenceTable(self):
		self.logger.error("Not yet implemented. do something!")

def convergenceFrontend(convergencetest, description=__doc__):
	"""
	Use this as your main function and pass an instance of the convergence test
	which shall be run after argument dispatching.
	"""
	logger = logging.getLogger("convergenceFrontend")

	frontend = ExaFrontend(program_description=description)
	reportAdapter = ConvergenceStarterReportingAdapter()

	frontend.add_module(convergencetest)
	frontend.add_module(reportAdapter)

	args = frontend.parse_args()
	processes = convergencetest.startFromArgs(frontend.args, frontend.parser)
	reportAdapter.dispatchProcesses(processes)

	if reportAdapter.reporting:
		logger.info("Finishing with convergencePassed=%s" % str(reportAdapter.convergencePassed))
		# todo: let this the convergence_table do.
		sys.exit(0 if reportAdapter.convergencePassed else -3)
	else:
		logger.info("All done.")


