#!/usr/bin/env python
#
# This script is part of the "Exa" Python Postprocessing tools.
# They are written in Python 2 with Python 3 compatibility in mind.
# Dependencies: VTK and numpy/scipy.
#
# This file is part of the ExaHyPE project.
# Copyright (c) 2016  http://exahype.eu
# All rights reserved.
#
# The project has received funding from the European Union's Horizon 
# 2020 research and innovation programme under grant agreement
# No 671698. For copyrights and licensing, please consult the webpage.
#
# Released under the BSD 3 Open Source License.
# For the full license text, see LICENSE.txt

"""
The ExaHelpers module collects utility code for the exareader.py
and the exaplayer.py scripts.
"""

from __future__ import print_function # Py3 compilance
import sys, argparse, logging
import numpy as  np # for vectorize_concatenate

logger = logging.getLogger("exahelpers")

# helper functions:
find_nearest = lambda a, a0: a.flat[abs(a-a0).argmin()]
unzip = lambda z: zip(*z) # inverse zip
middle = lambda a: a[int(len(a)/2)]
first = lambda l: l[0]
last = lambda l: l[-1]

class fileformat:
	"""
	A class representing a number of functions to read in or
	write out data. Formats are registered as functions via decorators.

	This class should have a more generic name as the concept is also
	quite generic and was already applied to grid formats for the vtk
	readers.
	"""
	def __init__(self, name):
		self.name=name
		self.formats={}
		self.description={}
	def register(self, command_name, desc=""):
		def decorator(func):
			self.formats[command_name] = func
			self.description[command_name] = desc
			return func
		return decorator
	def choices(self):
		return self.formats.keys()
	def get(self, chosenformat):
		"get() returns a function"
		if not chosenformat in self.formats:
			raise ValueError("Could not find format '%s' in %s" % (chosenformat, str(self)))
		logger.info("%s selection: %s" % (self.name, chosenformat))
		return self.formats[chosenformat]
	def desc(self, chosenformat):
		return self.description[chosenformat]
	def __str__(self):
		return 'Format choser class (%s) with %d formats registered: '+ \
			", ".join(["%s (%s)"%nameDesc for nameDesc in self.description.iteritems() ])


# convenient numpy method
def vectorize_concatenate(func):
	"""
	Decorator to vectorize file input readers like given in argio or vtkreader.
	"""
	def func_wrapper(*args, **kwargs):
		args_without_fname = list(args)
		fnames = args_without_fname.pop(0)
		if not isinstance(fnames ,basestring):
			# progress bar: works but not nice output
			#genlogwrap = lambda i: lambda text,*a,**kw: log("[%02d/%02d] "% (i,len(fname))+text,*a,**kw)
			
			# optimization to avoid np.vstack for large single input files (>3GB):
			if len(fnames) == 1:
				return func(fnames[0], *args_without_fname, **kwargs)

			outputs = [func(f, *args_without_fname, **kwargs) for i,f in enumerate(fnames)]
			# vstack does not really work (any more, for some reason) for my
			# recarrays.
			#return np.vstack(tuple(outputs))

			# instead, this does:
			import numpy.lib.recfunctions as nlr
			return nlr.merge_arrays(outputs, flatten=True)
		else:
			return func(*args, **kwargs)
	return func_wrapper

class ExaVerbosity:
	"""
	ExaVerbosity manages parsing and applying logging related arguments.
	"""
	def add_group(self, argparser):
		default_loglevel = logging.INFO
		group = argparser.add_argument_group('verbosity')
		group.add_argument('-q', '--quiet', action='store_const', dest='loglevel', default=default_loglevel,
		                   const=logging.WARNING, help="Quiet, do not give so many output")
		group.add_argument('-v', '--verbose', action='store_const', dest='loglevel', const=logging.DEBUG,
			           help='Be verbose, give more debugging output')

	def apply_args(self, args):
		logging.basicConfig(stream=sys.stderr, level=args.loglevel)

class ExaFrontend:
	"""
	ExaFrontend is just a tiny wrapper around argparse.ArgumentParser which allows neat module
	callbacks after parsing. And also adds the ExaVerbosity.
	"""
	def __init__(self, program_description="", epilog=""):
		programname = sys.argv[0]
		self.parser = argparse.ArgumentParser(description=program_description, epilog=epilog,
			formatter_class=argparse.RawDescriptionHelpFormatter)
		self.modules = []

	def add_module(self, argumentable):
		argumentable.add_group(self.parser)
		self.modules.append(argumentable)

	def parse_args(self, args=None):
		"""
		Calls argparser.parse_args. You can pass the `args` instead of `sys.argv`
		"""
		# add the logging stuff as the last group
		self.verbosity = ExaVerbosity()
		self.add_module(self.verbosity)

		self.args = self.parser.parse_args(args)
	
		# allow each module to do work after parsing
		for mod in self.modules:
			try:
				mod.apply_args(self.args)
			except AttributeError:
				pass

		return self.args

