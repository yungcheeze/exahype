#
# This is a small library very similar to the
#   Postprocessing/exahelpers.py
# library. It contains classes to glue together classes
# which are controlled by command line.
#

import logging, argparse, sys

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

	def apply_args(self, args, argparser):
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

	def add_module(self, argumentable, push_front=False):
		"""
		Allow argumentable to create a group or add arguments to the parser
		and then add it to our list of modules where to apply_args later on
		either at the beginning or end.
		"""
		argumentable.add_group(self.parser)
		if push_front:
			self.modules.insert(0, argumentable)
		else:
			self.modules.append(argumentable)

	def parse_args(self, args=None):
		"""
		Calls argparser.parse_args. You can pass the `args` instead of `sys.argv`
		"""
		# add the logging stuff as the last group displayed
		# but the first one being called apply_args on.
		self.verbosity = ExaVerbosity()
		self.add_module(self.verbosity, push_front=True)

		self.args = self.parser.parse_args(args)
	
		# allow each module to do work after parsing
		for mod in self.modules:
			try:
				mod.apply_args(self.args, self.parser)
			except AttributeError:
				pass

		return self.args
