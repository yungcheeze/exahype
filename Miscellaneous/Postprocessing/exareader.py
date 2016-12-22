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
exareader.py is both a standalone program as well a mini library to read VTK
ExaHyPE files as well as to convert them to various other file formats.

Given a typical set of files "solution-0.vtk, solution-1.vtk, .... solution-50.vtk, ..."
with each 200MB in size, containing vtkUnstructuredGrids. In order to deal
with them in Scientific python, the class `ExaReader` can convert them to numpy
structured arrays with the columns

   (x, y, z, time, Q0, Q1, Q2, ... QN)

where Q[n] is the n'th conserved variable in ExaHyPE. These arrays can be dumped
to files or be used interactively if you invoke this script with "python -i" or
just import it from the python command line.

Use the minimalistic interface like:

   from exareader import exareader
   data = exareader(glob('solution-*.vtk'))
   somethingnew = data['Q0'] + data['Q1'] / sqrt(data['x']**2 + data['y']**2)

You can also use the more powerful object oriented interface like:

   from exareader import ExaReader
   reader = ExaReader()
   data = reader.read_files(...)

This script is also part of the exa postprocessing tools, it plays nicely with
ExaPlayer and your own scripts.

Example invocations of the command line utility:

    ./exareader.py -o solution-0.txt.gz -c solution-0.vtk
    ./exareader.py solution-*.vtk | gzip > foobar.csv.gz

If you have files solution-1.vtk ... solution-15.vtk ... then you probably
want to ensure you glob with natural sorting, eg. with

    ./exareader.py -o data.csv.gz -c $(ls *vtk | sort -V)

ExaReader will never change the order of the input data. The output array is
basically just a vertical stacking of the data contained in the VTK files.
"""

from __future__ import print_function # Py3 compilance

# VTK libraries for Python; NumPy
from vtk import vtkUnstructuredGridReader
from vtk.util.numpy_support import vtk_to_numpy as npize # numpyize
import numpy as np

# Python built in batteries 2.7
import sys, logging, gzip
from collections import namedtuple

# our package
from exahelpers import fileformat, vectorize_concatenate, ExaFrontend

logger = logging.getLogger("exareader")

read_formats = fileformat('Input Fileformat')
write_formats = fileformat('Output Fileformat')
grid_formats = fileformat('Grid Format')

@grid_formats.register('vertices', desc='vtk::Cartesian::vertices format')
@vectorize_concatenate
def read_exahype_cartesian_vertices(fname):
	"""
	Transforms a single VTK files to a recorded numpy array.
	Currently only filenames given as string are supported, no python streams.
	For a 200MB input file, this function needs around 1-2sec.

	This function represents the way of reading in the "OLD" ExaHyPE output file format
	(VTK format) which was used up to ~ Aug 2016.
	Back then, the payload were stored in Point Data as individual scalar fields.
	With the new format, the VTK files store datas in Cell Data where all payload
	is stored in a single vector field.

	Probably this also interferes with the differences of
		the ExaHyPE vtk::Cartesian::cells::ascii plotter
		vs. the ExaHyPE vtk::Cartesian::vertices::ascii plotter
	
	This vtkreader is probably still suitable for the vertices plotter.
	It has not been tested, thought. [01. Oct 2016]	
	"""

	logger.info("Reading %s... " % fname)
	reader = vtkUnstructuredGridReader()
	reader.SetFileName(fname)
	reader.ReadAllScalarsOn()
	reader.Update()
	dataset = reader.GetOutput() # is a VtkUnstructuredGrid
	points = dataset.GetPoints()
	scalars = dataset.GetPointData()
	num_scalars = scalars.GetNumberOfArrays()
	classname = dataset.GetClassName()
	num_cells = dataset.GetNumberOfCells()
	num_points = dataset.GetNumberOfPoints()
	logger.info("File is a %s with %s cells, %s points, %d scalar fields" % (classname, num_cells, num_points, num_scalars))

	point_list = npize( points.GetData() )
	datasets = point_list.transpose()
	column_names = ['x','y','z']

	for i in range(0, num_scalars):
		scalar_field = npize( scalars.GetArray(i) )
		datasets = np.vstack((datasets,scalar_field))
		column_names.append(scalars.GetArrayName(i))

	logger.info("Exporting scalar fields "+str(column_names))

	ret = datasets.transpose()
	#numpy.float32
	structured_array_dtype = np.dtype([ (k, ret.dtype) for k in column_names ])
	ret.dtype = structured_array_dtype
	return ret

@grid_formats.register('cells', desc='vtk::Cartesian::cells format')
@vectorize_concatenate
def read_exahype_cartesian_cells(fname):
	"""
	Transforms a single VTK files to a recorded numpy array.
	This function is the logical successor of the old vtkreader(fname) function.
	While the former is probably still valid for the exahype::Cartesian::vertices
	format, it is definetly no more for the new one, see above.

	This rewrite took quite some time but it is also much faster then the original
	function (as it omits the loop over the individual fields, copying a huge amount
	of data around).
	"""
	logger.info("Reading according to exahype::cartesian::cells the file %s... " % fname)
	reader = vtkUnstructuredGridReader()
	reader.SetFileName(fname)
	reader.ReadAllScalarsOn()
	reader.Update()
	dataset = reader.GetOutput() # is a VtkUnstructuredGrid

	# we can get a lot of debugging information:
	logger.debug("All information about the VTK data structure: "+str(dataset))

	classname = dataset.GetClassName()
	num_cells = dataset.GetNumberOfCells()
	num_points = dataset.GetNumberOfPoints()

	logger.info("File is a %s with %s points and %s cells" % (classname, num_points, num_cells))

	# the new file format has no more point data, so this holds:
	pointdata = dataset.GetPointData()
	num_pointdata = pointdata.GetNumberOfArrays()
	assert num_pointdata == 0

	# access the fields stored in the cells
	celldata = dataset.GetCellData()
	numfields = celldata.GetNumberOfArrays()
	arrays = { celldata.GetArrayName(i): npize(celldata.GetArray(i)) for i in range(numfields) }
	for n,v in arrays.iteritems():
		logger.info("CellData array '%s' has the shape %s" % (n, str(v.shape)))

	# here we put the assumption that the interesting array has the name "Q"
	# and "time" holds the time information about each cell
	PayloadArrayName = "Q" 
	TimeArrayName = "time"
	Q = arrays[PayloadArrayName]
	Q_column_names = ["Q%d"%d for d in range(len(Q[0])) ]

	# resolve the coordinates of the cells
	# there are also fancy methods requiring a cellindex each:
	# 1) dataset.GetCellPoints
	#    to get a list of points which define the boundary of the cell.
	# 2) dataset.GetCellBounds(cellId, [xmin,xmax,ymin,ymax,zmin,zmax])
	#    to get the bounding box of the cell, as our cells are rectangular
	#    this is merely a shortcut of GetCellPoints.
	# 3) dataset.GetPointsCells
	#    the inverse of GetCellPoints: get a list of cells which are
	#    flanking a given point. The way to go to interpolate field values
	#    at the grid points.

	# we can get the locations of all points like that:
	points = dataset.GetPoints()
	point_locs = npize(points.GetData())

	# and we can get the cell connectivity information
	# it is in the list format  (n,id1,id2,...,idn, n,id1,id2,...,idn, ...)
	# with n is the number of points in the cell,
	#      id is a zero-offset index into an associated point list.
	cells = dataset.GetCells()
	conn = npize(cells.GetData())

	npoints = conn[0] # in 2D, we have vtkPixels, so this is 4
	spanwidth = npoints + 1 # periodicy in the list format, 1 == len([npoints])
	all_npoints = conn[::spanwidth]
	assert len(all_npoints) == num_cells # assert the list length
	assert np.all(all_npoints == npoints) # assert the list format

	# reshape and strip the leading '4' (npoints) column
	# conntable now holds for every cell the list of associated points
	# By inspecting point_locs[conntable], we can see for each cell
	# the list of coordinates of the points which define its boundaries.
	conntable = conn.reshape((num_cells, spanwidth))[:, 1:]
	
	# we assign each cell a coordinate by just taking the *minimum* coordinate
	# of each point belonging to the cell. Now we have the well-known
	# (x,y,z) array shaped (num_cells, 3).
	cell_locs = np.amin(point_locs[ conntable ], axis=1)
	cell_locs_column_names = list("xyz")

	# prepare time as list of columns, suitable for merging
	time = arrays[TimeArrayName][:, np.newaxis]
	time_column_names = ["time"]

	# now merge all the columns to the good old structured dataset
	merged_bigtable = np.hstack((time, cell_locs, Q))
	column_names = time_column_names + cell_locs_column_names + Q_column_names

	def structured_conversion_magic(nparray, colnames):
		npt = nparray.T
		npt.dtype = np.dtype([ (k, nparray.dtype) for k in colnames ])
		return npt[0] # why the heck

	return structured_conversion_magic(merged_bigtable, column_names)

@write_formats.register('np', desc='Binary numpy file')
def output_npy(npdata, outputfname, **args):
	np.save(outputfname, npdata)

@write_formats.register('csv', desc='Comma-Seperated-Values')
def output_ascii(npdata, outputfname, **args):
	logger.info("This may take a while. Have a coffee. Or watch your output growing.")
	np.savetxt(
		fname=outputfname,
		X=npdata,
		fmt=args['nformat'], # args type (dict/Namespace) not really clear here
		header=','.join(['"%s"' % s for s in npdata.dtype.names]),
		comments=''
	)

@read_formats.register('vtk', desc='Vizualisation Toolkit file')
def input_vtk(filenames, gridtype='cells', **args):
	logger.info("Invoking VTK interface with %d filenames" % len(filenames))
	vtkreader = grid_formats.get(gridtype)
	return vtkreader(filenames)

@read_formats.register('np', desc='Binary numpy file')
@vectorize_concatenate
def input_numpy(filename, **args):
	logger.info("Loading numpy file %s" % filename)
	return np.load(filename)

@read_formats.register('auto', desc='Auto detect')
def input_autodetect(filenames, **args):
	"Detect the format of inputfiles based on the filenames"
	
	# I know os.path.splitext would be cleaner, but this approach works, too
	detectors = {
		"vtk": lambda fname: ".vtk" in fname,
		"np":  lambda fname: ".np" in fname
	}
	for fileformat, detector in detectors.iteritems():
		detectormap = map(detector, filenames)
		if all(detectormap):
			fundamentalreader = read_formats.get(fileformat, 'autodetected')
			return fundamentalreader(filenames, **args)
		if any(detectormap):
			logger.info("Detected %d files of %s format" % (sum(detectormap), fileformat))
			
	raise LookupError("Could not determine file format inputfiles based on their filenames. Please specify.")


class ExaReader:
	"""
	The reader class for managing read access to files
	"""
	def __init__(self, files_as_main=True, inputfiles_required=True):
		self.files_as_main = files_as_main
		self.inputfiles_required = inputfiles_required
		
		self.default_inputformat = 'auto'

	def add_group(self, argparser):
		group = argparser.add_argument_group('input')

		files_help = 'The file(s) to read in. Multiple VTK input files are simply attached in output. Also numpy files supported'
		if self.files_as_main:
			group.add_argument('inputfiles', metavar='solution-0.vtk', nargs='+', help=files_help)
		else:
			group.add_argument('-r', '--inputfiles', metavar='solution-0.vtk', nargs='+', required=self.inputfiles_required, help=files_help) # action='append',
		group.add_argument('-i', '--inform', dest='inform', choices=read_formats.choices(), type=str, default=self.default_inputformat,
		               help='File format of the input files, VTK takes long, numpy is rather quick.')
		group.add_argument('-g', '--gridtype', dest='gridtype', choices=grid_formats.choices(), type=str, default='cells',
		               help='Plotter format used during ExaHyPE run (staggered grid?)')

	def remove_prog_in_inputfiles_list(self, argparser, args):
		"""
		A workaround to remove the program name in the args.inputfiles, if present.
		It's not exact but kind of heuristic.
		"""
		logger.info("Removing program name '%s' from inputfiles list" % argparser.prog)
		programname = argparser.prog
		argparser.inputfiles = [a for a in args.inputfiles if not programname in a]
		return argparser

	def read_files_as_requested(self, args):
		"""
		Convenience call to convert argparse.Namespace objects to method call parameters.
		"""
		return self.read_files(**vars(args))

	def read_files(self, inform="np", inputfiles=[], **readerkwargs):
		"""
		Return the numpy array which is constructed according to arguments.
		"""
		if not inputfiles or not len(inputfiles):
			raise ValueError("Please provide at least one input file.")
		reader = read_formats.get(inform)
		logger.info("Reading input as %s from the following files:" % read_formats.desc(inform))
		for i,f in enumerate(inputfiles): logger.info(" %d. %s" % (i+1,f))
		data = reader(inputfiles, **readerkwargs)
		return data

def exareader(fnames, quiet=False):
	"""
	This is back the pythonic function interface how to read in data
	"""
	if not quiet:
		# enable logging, as this is an entrypoint into the module
		logging.basicConfig(stream=sys.stderr, level=args.INFO)
	reader = ExaReader()
	return reader.read_files(inputfiles=fnames, inform="vtk")

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

		logger.info("Printing output to %s" % str(outputfh.name))
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

def exatractor(args=None):
	"""
	The Command Line Interface, can also be used as a python function. While args
	are the command line arguments for the parser, typically you can give it just
	the names of your input files as this is the default of ExaReader
	(cf. ExaReader.files_as_main). That is, an example usage is

	>>> exareader(glob('*.vtk'))
	"""
	frontend = ExaFrontend(epilog=__doc__,
		program_description='ExaHyPE simulation data converter, ie. VTK to plain ASCII converter')

	reader = ExaReader(files_as_main=True)
	writer = ExaWriter()

	frontend.add_module(reader)
	frontend.add_module(writer)
	args = frontend.parse_args(args)

	logger.info("Welcome to the ExaReader/ExaWriter CLI")
	data = reader.read_files_as_requested(args)
	logger.info("Have read a %s-shaped numpy array", data.shape)
	writer.write_output_as_requested(data, args)
	logger.info("Finished")

if __name__ == "__main__":
	exatractor()


