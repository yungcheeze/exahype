#!/usr/bin/python
#
# This script is part of ExaHyPE. Probably.
# This script is Python 2.
#
# This is a small helper script to convert VTK ExaHype files to plain ASCII. These
# are typically a set of files "solution-0.vtk, solution-1.vtk, ... solution-50.vtk, ..."
# with each 200MB in size, containing vtkUnstructuredGrids. 
# This script collects them and produces a single output file which may be ASCII
# for convenience.
#
# This script is not performant at all, it collects all input in numpy arrays before
# writing any output. In a typical scenario with 40 output files, 3GB in size, the
# script needs around 8 mins and 2GB system memory before writing a 200MB compressed
# output file.
#
# -- Public Domain, May, 10. 2016, Sven K.

from __future__ import print_function

from vtk import *
from vtk.util import numpy_support
import numpy as np

# Python built in batteries
import argparse, sys # Python 2.7
from collections import namedtuple # Python 2.4
import gzip

HeaderedTable = namedtuple('HeaderedTable', ['column_names', 'data'])

verbose = False
def log(text, newline=True, force=True):
	if not force and not verbose: return
	print(text, file=sys.stderr, end="\n" if newline else '')
	sys.stderr.flush()

def vtk2data(fname, verbose=False):
	"""
	Transforms a single VTK file to a numpy array (+ column names)
	"""
	log("Reading %s... " % fname, newline=False)
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
	log("is a %s with %s cells, %s points, %d scalar fields" % (classname, num_cells, num_points, num_scalars))

	point_list = numpy_support.vtk_to_numpy( points.GetData() )
	datasets = point_list.transpose()
	column_names = ['x','y','z']

	for i in range(0, num_scalars-1):
		scalar_field = numpy_support.vtk_to_numpy( scalars.GetArray(i) )
		datasets = np.vstack((datasets,scalar_field))
		column_names.append(scalars.GetArrayName(i))

	return HeaderedTable(column_names, datasets.transpose())

def vtks2data(fnames, verbose=False):
	"""
	Transforms multiple VTK files to a single numpy array (+ column names)
	"""
	outputs = [vtk2data(f, verbose) for f in fnames]
	datasets = np.vstack(tuple([ht.data for ht in outputs]))
	column_names = outputs[0].column_names
	return HeaderedTable(column_names, datasets)

def ht2output(headeredtable, outputfname, fmt):
	np.savetxt(
		fname=outputfname,
		X=headeredtable.data,
		fmt=fmt,
		header=','.join(['"%s"' % s for s in headeredtable.column_names]),
		comments=''
	)


if __name__ == "__main__":
	programname = sys.argv[0]
	parser = argparse.ArgumentParser(description='ExaHyPE VTK to plain ASCII converter',
		epilog=("Example invocations:",
			"   %s -o solution-0.txt.gz -c solution-0.vtk" % programname,
			"   %s solution-*.vtk | gzip > foobar.csv.gz" % programname,
			"if you have files solution-1.vtk ... solution-15.vtk ... than you probably",
			"want natural sorting, eg. with",
			"   %s -o data.csv.gz -c $(ls *vtk | sort -V)" % programname,
			"%s will not never change the order of the input data." % programname
		))

	parser.add_argument('fnames', metavar='solution-0.vtk', nargs='+',
	                   help='The VTK file(s) to convert. Multiple input files are simply attached in output.')
	parser.add_argument('-c', '--compress', dest='compress', action='store_true', default=False,
	                   help='Compress the output, creates a gzip file')
	parser.add_argument('-q', '--quiet', dest='quiet', action='store_true', default=False,
	                   help='Quiet, no informal output on stderr')
	parser.add_argument('-o', '--output', dest='output', default=False,
	                   help='Output filename (default: stdout)')
	parser.add_argument('-f', '--format', dest='format', default='%.5e',
			   help='Number format string, cf. numpys savetxt() documentation')

	args = parser.parse_args()
	verbose = not args.quiet
	use_stdout = not args.output

	opener = gzip.open if args.compress else open
	if args.compress and use_stdout:
		log("Gzip output for stdout currently not supported. Just use pipes: %s ... | gzip" % programname, force=True)
		sys.exit(1)
	outputfh = sys.stdout if not args.output else opener(args.output, 'w')

	headeredtable = vtks2data(args.fnames, verbose)
	log("Printing output to %s" % str(outputfh))
	try:
		ht2output(headeredtable, outputfh, args.format)
	except IOError:
		if use_stdout:
			# eg. when using this command | head
			# stdout is closed, no point in continuing
			# Attempt to close them explicitly to prevent cleanup problems:
			try:
				sys.stdout.close()
			except IOError:
				pass
			try:
				sys.stderr.close()
			except IOError:
				pass
	log("Finished")

