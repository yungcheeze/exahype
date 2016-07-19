# python2 lib for argparse IO helpers

from __future__ import print_function
import numpy as np

# same module
import vtkreader
from . import vectorize_concatenate, log

# module variables
mylog = log

class fileformat:
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
		# returns a function
		return self.formats[chosenformat]
	def desc(self, chosenformat):
		return self.description[chosenformat]

informats = fileformat('input')
outformats = fileformat('output')

@outformats.register('np', desc='Binary numpy file')
def output_npy(npdata, outputfname, args=None):
	np.save(outputfname, npdata)

@outformats.register('csv', desc='Comma-Seperated-Values')
def output_ascii(npdata, outputfname, args):
	mylog("This may take a while. Have a coffee. Or watch your output growing.")
	# TODO: I inserted a bug here. Doesn't work any more with ASCII here. Dunno why.
	np.savetxt(
		fname=outputfname,
		X=npdata,
		fmt=args.nformat,
		header=','.join(['"%s"' % s for s in npdata.dtype.names]),
		comments=''
	)

@informats.register('vtk', desc='Vizualisation Toolkit file')
def input_vtk(filenames, args=None):
	logger = logger_for(args) if args else log
	return vtkreader.vtkreader(filenames, logger)

@informats.register('np', desc='Binary numpy file')
@vectorize_concatenate
def input_numpy(filename, args=None):
	log("Loading numpy file %s" % filename)
	return np.load(filename)


def add_io_group(argparser):
	group = argparser.add_argument_group('input/output')
	
	group.add_argument('inputfiles', metavar='solution-0.vtk', nargs='+',
                           help='The file(s) to read in. Multiple VTK input files are simply attached in output. Also numpy files supported')
	group.add_argument('-i', '--inform', dest='inform', choices=informats.choices(), type=str, default='vtk',
		           help='File format of the input files, VTK takes long, numpy is rather quick.')
	group.add_argument('-c', '--compress', dest='compress', action='store_true', default=False,
	                   help='Compress the output, creates a gzip file. Applies only if output is not stdout.')
	group.add_argument('-q', '--quiet', dest='quiet', action='store_true', default=False,
	                   help='Quiet, no informal output on stderr')
	group.add_argument('-f', '--outfile', dest='outfile', default=False,
	                   help='Output filename (default: stdout)')
	group.add_argument('-o', '--outform', dest='outform', choices=outformats.choices(), type=str, default='csv',
			   help='Output file formats. CSV (ASCII) takes long, binary formats are rather quick.')
	group.add_argument('-n', '--numberformat', dest='nformat', default='%.5e',
			   help='Number format string, cf. numpys savetxt() documentation. Applies only for CSV.')


def logger_for(parsed_args):
	def emptylogger(text, newline=True, force=True): pass
	logger = emptylogger if parsed_args.quiet else vtkreader.log
	mylog = logger
	return logger

def get_input(args):
	"""
	Return the numpy array which is constructed according to the input formats.
	"""
	log = logger_for(args)
	reader = informats.get(args.inform)
	log("Reading input as %s from %s" % (informats.desc(args.inform), args.inputfiles))
	data = reader(args.inputfiles, args)
	return data

def write_output(data, args):
	log = logger_for(args)
	use_stdout = not args.outfile
	opener = gzip.open if args.compress else open
	# this parameter check is basically too late.
	if args.compress and use_stdout:
		raise AttributeError("Gzip output for stdout currently not supported. Just use pipes: %s ... | gzip" % programname)
	outputfh = sys.stdout if not args.outfile else opener(args.outfile, 'w')

	log("Printing output to %s" % str(outputfh))
	try:
		writer = outformats.get(args.outform)
		writer(data, outputfh, args)
	except IOError:
		if use_stdout:
			# eg. when calling python, calling this method and piping on the shell to "./python ... | head",
			# then stdout is closed, no point in continuing
			# Attempt to close them explicitly to prevent cleanup problems:
			try:
				sys.stdout.close()
			except IOError:
				pass
			try:
				sys.stderr.close()
			except IOError:
				pass
