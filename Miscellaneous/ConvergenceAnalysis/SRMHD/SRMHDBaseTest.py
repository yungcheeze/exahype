#!/usr/bin/env python
#
# Run convergence tests.
#
# This is basically python version agnostic and should run with both python2
# and python3.
#
# ExaHyPE 2016, SvenK
#

"""
Run convergence tests: The SRMHD tests.
Will setup environment variables to execute a templated specification file
If called without arguments, it will run a series of runs in parallel.

Sample usages from the command line:

$ ./AlfenWaveTest.py run -p 2 -m 0.1
$ ./AlfenWaveTest.py run-report --all
"""

import sys, logging
logger = logging.getLogger("runSRMHDTest")
sys.path.append("../libconvergence")

from convergence_test import PolyorderTest
from convergence_frontend import ConvergenceFrontend

class SRMHDTest(PolyorderTest):
	"""
	The Ideal Special Relativistic MHD test is a PolyorderTest is a convergence
	test on a uniform grid where both the grid spacing as well as the polynomial
	order in ADERDG is changed in order to see convergence accross both degrees
	of freedom (with reasonable SRMHD defaults).
	"""
	def __init__(self, name="UnnamedSRMHDTest"):
		PolyorderTest.__init__(self, name)

		self.settings['ExaBinary'] = "../../../ApplicationExamples/SRMHD/ExaHyPE-MHDSolver"
		self.settings['ExaSpecfile'] = "MHD_GenericConvergence.exahype"
		# template to set up a queueing system
		self.settings['QRUNTPL'] = "" #"srun -n1 --partition=x-men --time=29:00:00 --mem=0 --job-name=p{ExapOrder}-m{ExaMeshSize}-AlfenWave"

		# parameters for setting up the specfile
		self.settings['ExaEndTime'] = 10.0
		# parameters deciding how frequently output is made. As a first criterion,
		# 1 output dump with the highest resolution is 250MB.
		self.settings['ExaConvOutputRepeat'] = 0.1
		self.settings['ExaVtkOutputRepeat'] = 0.2
		# single threaded in the moment.
		self.settings['ExaTbbCores'] = 1

		# this is useful if you compiled with assertions
		self.settings['EXAHYPE_SKIP_TESTS'] = True


	def setInitialData(self, name):
		# options are: {"alfenwave", "blast", "orsagtang", "rotor"}
		self.settings['EXAHYPE_INITIALDATA'] = name

	def setExtend(self, width, height):
		self.polyorders = range(2,10)
		self.width = width
		self.height = height
		self.depths = range(1,6)

		self.computeCombinations()
		until = self.until # shorthand

		self.adaptiveRunRange = {
			2: until(243), 3: until(243), 4: until(243),
			5: until(243), 6: until(81), 7: until(81),
			8: until(27), 9: until(27)
		}

		self.settings['ExaWidth'] = self.width
		self.settings['ExaHeight'] = self.height



