#!/usr/bin/env python
#
# This script is part of ExaHyPE. Probably.
# This script is Python 2.
#
# This python script is rather a mini library to read VTK ExaHyPE files,
# given as atypical set of files "solution-0.vtk, solution-1.vtk, .... solution-50.vtk, ..."
# with each 200MB in size, containing vtkUnstructuredGrids, to numpy structured
# numpy arrays with the columns
#
# (x, y, z, time, Q0, Q1, Q2, ... QN)
#
# where Q[n] is the n'th conserved variable in ExaHyPE.
#
# Use it like:
#
#   from vtkreader import vtkreader
#   data = vtkreader(glob('solution-*.vtk'))
#   somethingnew = data['Q0'] + data['Q1'] / sqrt(data['x']**2 + data['y']**2)
#
# This script is also part of the exa postprocessing tools.
# 
#
# -- Public Domain, 2016, Sven K.

from __future__ import print_function # Py3 compilance

# VTK libraries for Python; NumPy
from vtk import vtkUnstructuredGridReader
from vtk.util import numpy_support
import numpy as np

# Python built in batteries 2.7
import argparse, sys, logging, gzip
from collections import namedtuple

logger = logging.getLogger("exareader")

# convenient python2 methods
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
			
			# optimization to avoid np.vstack for large input files (>3GB):
			if len(fnames) == 1:
				return func(fnames[0], *args_without_fname, **kwargs)

			outputs = [func(f, *args_without_fname, **kwargs) for i,f in enumerate(fnames)]
			return np.vstack(tuple(outputs))
		else:
			return func(*args, **kwargs)
	return func_wrapper

@vectorize_concatenate
def vtkreader(fname):
	"""
	Transforms a single VTK files to a recorded numpy array.
	Currently only filenames given as string are supported, no python streams.
	For a 200MB input file, this function needs around 1-2sec.
	"""

	logger.info("Reading %s... " % fname)
	reader = vtkUnstructuredGridReader()
	reader.SetFileName(fname)
	reader.ReadAllScalarsOn()
	reader.Update()
	dataset = reader.GetOutput()
	points = dataset.GetPoints()
	scalars = dataset.GetPointData()
	num_scalars = scalars.GetNumberOfArrays()
	classname = dataset.GetClassName()
	num_cells = dataset.GetNumberOfCells()
	num_points = dataset.GetNumberOfPoints()
	logger.info("File is a %s with %s cells, %s points, %d scalar fields" % (classname, num_cells, num_points, num_scalars))

	point_list = numpy_support.vtk_to_numpy( points.GetData() )
	datasets = point_list.transpose()
	column_names = ['x','y','z']

	for i in range(0, num_scalars):
		scalar_field = numpy_support.vtk_to_numpy( scalars.GetArray(i) )
		datasets = np.vstack((datasets,scalar_field))
		column_names.append(scalars.GetArrayName(i))

	logger.info("Exporting scalar fields "+str(column_names))

	ret = datasets.transpose()
	#numpy.float32
	structured_array_dtype = np.dtype([ (k, ret.dtype) for k in column_names ])
	ret.dtype = structured_array_dtype
	return ret

class fileformat:
	"""
	A class representing a number of functions to read in or
	write out data. Formats are registered as functions via decorators.
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
		# returns a function
		logger.info("Fileformat selection: %s" % chosenformat)
		return self.formats[chosenformat]
	def desc(self, chosenformat):
		return self.description[chosenformat]

read_formats = fileformat('input')
write_formats = fileformat('output')

@write_formats.register('np', desc='Binary numpy file')
def output_npy(npdata, outputfname, **args):
	np.save(outputfname, npdata)

@write_formats.register('csv', desc='Comma-Seperated-Values')
def output_ascii(npdata, outputfname, **args):
	logger.info("This may take a while. Have a coffee. Or watch your output growing.")
	# TODO: I inserted a bug here. Doesn't work any more with ASCII here. Dunno why.
	np.savetxt(
		fname=outputfname,
		X=npdata,
		fmt=args['nformat'], # args type (dict/Namespace) not really clear here
		header=','.join(['"%s"' % s for s in npdata.dtype.names]),
		comments=''
	)

@read_formats.register('vtk', desc='Vizualisation Toolkit file')
def input_vtk(filenames, **args):
	logger.info("Invoking VTK interface with %d filenames" % len(filenames))
	return vtkreader(filenames)

@read_formats.register('np', desc='Binary numpy file')
@vectorize_concatenate
def input_numpy(filename, **args):
	logger.info("Loading numpy file %s" % filename)
	return np.load(filename)


class ExaReader:
	"""
	The reader class for managing read access to files
	"""
	def __init__(self, files_as_main=True):
		self.files_as_main = files_as_main

	def add_group(self, argparser):
		group = argparser.add_argument_group('input')

		files_help = 'The file(s) to read in. Multiple VTK input files are simply attached in output. Also numpy files supported'
		if self.files_as_main:
			group.add_argument('inputfiles', metavar='solution-0.vtk', nargs='+', help=files_help)
		else:
			group.add_argument('-r', '--inputfiles', metavar='solution-0.vtk', action='append', help=files_help)
		group.add_argument('-i', '--inform', dest='inform', choices=read_formats.choices(), type=str, default='vtk',
				   help='File format of the input files, VTK takes long, numpy is rather quick.')

	def read_files_as_requested(self, args):
		"""Convenience call to convert argparse.Namespace objects to method call parameters"""
		return self.read_files(**vars(args))

	def read_files(self, inform="np", inputfiles=[], **readerargs):
		"""
		Return the numpy array which is constructed according to arguments.
		`args` is a 
		"""
		reader = read_formats.get(inform)
		logger.info("Reading input as %s from the following files:" % read_formats.desc(inform))
		for i,f in enumerate(inputfiles): logger.info(" %d. %s" % (i+1,f))
		data = reader(inputfiles, **readerargs)
		return data


class ExaWriter:
	"""
	The writer class doing the same as reader just for writing.
	"""
	def add_group(self, argparser):
		group = argparser.add_argument_group('output')
		group.add_argument('-c', '--compress', dest='compress', action='store_true', default=False,
			           help='Compress the output, creates a gzip file. Applies only if output is not stdout.')
		group.add_argument('-f', '--outfile', dest='outfile', default=False,
			           help='Output filename (default: stdout)')
		group.add_argument('-o', '--outform', dest='outform', choices=write_formats.choices(), type=str, default='csv',
				   help='Output file formats. CSV (ASCII) takes long, binary formats are rather quick.')
		group.add_argument('-n', '--numberformat', dest='nformat', default='%.5e',
				   help='Number format string, cf. numpys savetxt() documentation. Applies only for CSV.')

	def write_output_as_requested(self, data, args):
		"""
		Convenience call to convert argparse.Namespace objects to method call parameters.
		Use this function if you have collected your arguments by argparser.
		Else use the pythonic method call `write_output` directly.
		"""
		return self.write_output(data, **vars(args))

	def write_output(self, data, outfile, outform="csv", compress=False, **writerargs):
		"""
		Writes out `data` to `outfile` using the global `write_formats` class where
		writers are registered.

		Parameters are the same as described by the `add_group` method.
		Furthermore, for the parameters without default arguments:
		data: Data to write in numpy format
		outfile: filename where to write to. If 'False' or 'None', will write to stdout.		
		"""
		use_stdout = not outfile
		opener = gzip.open if compress else open
		# this parameter check is basically too late.
		if compress and use_stdout:
			raise AttributeError("Gzip output for stdout currently not supported. Just use pipes: %s ... | gzip" % programname)
		outputfh = sys.stdout if not outfile else opener(outfile, 'w')

		logger.info("Printing output to %s" % str(outputfh))
		try:
			writer = write_formats.get(outform)
			writer(data, outputfh, **writerargs)
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


def vtkextractor():
	programname = sys.argv[0]
	logging.basicConfig(stream=sys.stderr, level=logging.DEBUG)
	parser = argparse.ArgumentParser(description='ExaHyPE simulation data converter, ie. VTK to plain ASCII converter',
		epilog="\n".join(["Example invocations:",
			"   %s -o solution-0.txt.gz -c solution-0.vtk" % programname,
			"   %s solution-*.vtk | gzip > foobar.csv.gz" % programname,
			"if you have files solution-1.vtk ... solution-15.vtk ... than you probably",
			"want natural sorting, eg. with",
			"   %s -o data.csv.gz -c $(ls *vtk | sort -V)" % programname,
			"%s will not never change the order of the input data." % programname
		]),
		formatter_class=argparse.RawDescriptionHelpFormatter)

	logger.info("Welcome to the ExaReader/ExaWriter CLI (called %s)" % programname)
	reader = ExaReader()
	writer = ExaWriter()

	reader.add_group(parser)
	writer.add_group(parser)
	args = parser.parse_args()

	data = reader.read_files_as_requested(args)
	logger.info("Have read a %s-shaped numpy array", data.shape)
	writer.write_output_as_requested(data, args)
	logger.info("Finished")


if __name__ == "__main__":
	vtkextractor()


