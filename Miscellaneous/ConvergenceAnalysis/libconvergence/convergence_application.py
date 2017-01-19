#!/usr/bin/env python
# Run convergence tests.
#
# This is basically python version agnostic and should run with both python2
# and python3.
#
# ExaHYPE 2016, SvenK

"""
The Convergence application package contains a main class to manage a
convergence test with a single frontend. The usage is quite simple:

> test = MyConvergenceTest("ShuVortexWhatever")
> test....()
> app = ConvergenceApplication(test, description="The program description")
> 
> if __name__ == "__main__":
> 	app.parse_args()

or alternatively, also from command line, like:

> app.parse_args(["report", "-a"])
"""

from numpy import arange, array
import subprocess, os, sys, time, logging # batteries

# libconvergence, same directory:
from convergence_arghelpers import ExaFrontend
from convergence_table import ConvergenceReporter, exitConvergenceStatus
from convergence_helpers import shell, getenv, pidlist, runBinary, cleandoc, MethodActions

class ConvergenceApplication:
	"""
	The frontend allows to link together the ConvergenceTest with the
	ConvergenceReporter. You can choose what the application should do by	
	choosing one of the registered actions.
	"""

	actions = MethodActions()

	def __init__(self, convergencetest, description=""):
		"""
		Use this as your main function and pass an instance of the convergence test
		which shall be run after argument dispatching.
		"""
		self.logger = logging.getLogger("ConvergenceFrontend")
		self.reporter = ConvergenceReporter()
		self.choices = self.actions.list()
		self.convergencetest = convergencetest

		self.frontend = ExaFrontend(program_description=description)
		self.frontend.add_module(self.convergencetest)
		self.frontend.add_module(self)
		self.frontend.add_module(self.reporter)

	def parse_args(self, args=None):
		"""
		This is the main method to start parsing and running.
		It basically passes the argument dispatching to ExaFrontend
		and then calls the appropriate action (run, report, ...).
		"""
		args = self.frontend.parse_args(args)
		self.actions.call(args.action, self)

	def add_group(self, parser):
		class_description = cleandoc(self)
		class_description += "\nThere are different modes how to run:\n\n"
		class_description += "\n".join([ "%10s: %s" % (action,doc) for (action,doc) in self.actions.docs.iteritems() ])

		group = parser.add_argument_group(title="ConvergenceFrontend", description=class_description)
		group.add_argument('action', default=None, choices=self.choices, help="Do one of the program actions")

	def waitForProcesses(self, processes):
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
		return exitcode

	@actions.add("run")
	def startConvergenceTest(self):
		"Start only the convergence test."
		self.logger.info("Starting convergence test")
		processes = self.convergencetest.startFromArgs(self.frontend.args, self.frontend.parser)
		return processes

	@actions.add("report")
	def startConvergenceReporter(self):
		"Start only the convergence report."
		self.logger.info("Starting convergence reporting")
		convergencePassed = self.reporter.start()
		self.logger.info("Finishing with convergencePassed=%s" % str(convergencePassed))
		exitConvergenceStatus(convergencePassed)

	@actions.add("run-wait")
	def startConvergenceTestAndWait(self):
		"Convergence test and wait until finished."
		processes = self.startConvergenceTest()
		exitcode = self.waitForProcesses(processes)
		return exitcode

	@actions.add("run-report")
	def startConvergenceTestAndWaitAndReporter(self):
		"Convergence test, wait until finished and do report afterwards."
		exitcode = self.startConvergenceTestAndWait()
		# important:
		# if exitcode == 0 or so!
		self.startConvergenceReporter()
		
