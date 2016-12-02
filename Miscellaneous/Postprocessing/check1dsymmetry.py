#!/usr/bin/env python
#

"""
This script checks in 2D and 3D simulations if the data is actually only 1D.
To do so, it uses the exareader infrastructure.
Typical simulations we run in exahype are plane waves, shocktubes, etc. which
are only one dimensional along one axis. We check this here.

Actually the clean way to do this would be an SVD / PCA to detect 1-dimensionality
in an arbitary direction.

This program does not plot anything, only CLI.
"""


from __future__ import print_function
import sys, argparse, logging

import numpy as np

from exahelpers import find_nearest, unzip, middle, first, last
from exahelpers import ExaFrontend
from exareader import ExaReader

logger = logging.getLogger("check1dsymmetry")

globalc = None

timecoord = "time"
spacecoords = ["x", "y", "z"]
allcoords = ["time", "x", "y", "z"]

class Check1DSymmetry:
	def __init__(self, sol, axis):
		self.sol = sol
		self.axis = axis
		self.uni = lambda col: np.unique(sol[col])

		self.times = self.uni(timecoord)
		self.x, self.y, self.z = map(self.uni, spacecoords)
		self.spacedim = 3 if len(self.z) > 1 else 2

		# Convert coordinates to continous array in order to search positions of
		# elements in solution.
		self.coords = sol[allcoords]
		# coords is now a recarray with (timesteps, number of points) rows
		self.coordtype = self.coords.dtype[0] # typically '<f4'
		self.nTimesteps = len(self.times)
		self.nPoints = self.coords.shape[1] / len(allcoords)
		assert self.nTimesteps == len(self.coords)
		# coords is now a continous array with the shape (timesteps, npoints, 4)
		self.coords = self.coords.view(self.coordtype).reshape(self.nTimesteps, self.nPoints, len(allcoords))

		""" AAARGH!:
In [197]: self.coords.size
Out[197]: 729000

In [198]: self.nTimesteps * self.nPoints * 4
Out[198]: 728960
		"""

	def passes(self):
		"Performs the actual check and returns true if success, false otherwise"
		# first: Setup an arbitrary array holding x,y,z data next to each other
		# instead of addressed like (t,x,y,z,Q0,Q1,Q2,...).
		# for each time.
		 
		

	@staticmethod
	def add_group(argparser):
		group = argparser.add_argument_group('check1dsymmetry')

		group.add_argument('-x', '--axis', choices=spacecoords, default="x",
			           help="On which axis to expect varying data (all other axis shall be equal)")
		group.add_argument('-t', '--tolerance', type=float, default=1e-10,
			           help="Deviation allowed")

	@staticmethod
	def main(args=None):
		frontend = ExaFrontend(epilog=__doc__,
			program_description='ExaHyPE 1d symmetry checker')

		reader = ExaReader(files_as_main=False)
		frontend.add_module(reader)
		frontend.add_module(Check1DSymmetry)

		args = frontend.parse_args(args)

		logger.info("Checking 1D symmetry")
		sol = reader.read_files_as_requested(args)
		checker = Check1DSymmetry(sol, args.axis)
		global globalc
		globalc = checker
		passed = checker.passes()
		if passed:
			logger.info("Check passed successfully")
			sys.exit(0)
		else:
			logger.error("Check failed")
			sys.exit(-1)



if __name__ == "__main__":
	Check1DSymmetry.main()

