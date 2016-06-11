#!/usr/bin/python
#
# This script is part of ExaHyPE. Probably.
# This script is Python 2.
#
# A quick way of "querying" data, or actually slicing data. Can be used both
# as a library or front end program.
#
# Currently only 1D slices are supported.
#

from __future__ import print_function # Py3 compilance
from vtkreader import vtkreader, log
import numpy as np
import sys, argparse

def slice_info(xyz, axis):
	search = dict(zip("xyz", list(xyz)))
	if not axis in search.keys():
		raise ValueError("Axis must be one out of %s" % str(search.keys()))
	del search[axis]
	return search

def slice_axis(data, xyz, axis):
	"""
	A very simple array slicer for numpy, based on coordinate values. Usage example:
		solution = vtkreader(glob("solutions-*.vtk"))
		xdata = slice_axis(solution, [ 0, 0.5, 0.5 ], 'x')
		xdata['x'] will be an increasing function
		xdata['y'] will be zero everywhere, same for xdata['z'].
		xdata.shape will less orders of magnitude smaller than solution.shape
	"""
	search = slice_info(xyz, axis)
	axis_filter = np.logical_and.reduce([data[ax]==v for ax,v in search.iteritems()])
	return data[axis_filter]


if __name__ == '__main__':
	programname = sys.argv[0]
	parser = argparse.ArgumentParser(description='ExaHyPE VTK data slicer to ASCII',
		epilog="blablabla")

	parser.add_argument('fnames', metavar='solution-0.vtk', nargs='+',
			   help='The VTK file(s) to convert. Multiple input files are simply attached in output.')
	parser.add_argument('-o', '--output', dest='output', default=False,
			   help='Output filename (default: stdout)')
	parser.add_argument('-n', '--numberformat', dest='nformat', default='%.5e',
			   help='Number format string, cf. numpys savetxt() documentation. Applies only for CSV.')

	queries = parser.add_argument_group('queries', 'Slicing description')
	#queries.add_argument('-X', '--xc', dest='x', default=None, help='Data point for X')
	#queries.add_argument('-Y', '--yc', dest='y', default=None, help='Data point for Y')
	#queries.add_argument('-Z', '--zc', dest='z', default=None, help='Data point for Z')
	queries.add_argument('-x', '--xp', dest='x', default=None, help='Grid (physical) point for X', type=float)
	queries.add_argument('-y', '--yp', dest='y', default=None, help='Grid (physical) point for Y', type=float)
	queries.add_argument('-z', '--zp', dest='z', default=None, help='Grid (physical) point for Z', type=float)

	queries.add_argument('-a', '--axis', dest='axis', choices="xyz", default='x', help="Which axis to slice.")

	args = parser.parse_args()
	use_stdout = not args.output

	outputfh = sys.stdout if not args.output else open(args.output, 'w')

	fulldata = vtkreader(args.fnames)

	xyz = [args.x, args.y, args.z]
	axis = args.axis
	description = "vtkslicer.py slices in %s-axis for %s" % (axis, slice_info(xyz, axis))
	log(description)
	sliced = slice_axis(fulldata, xyz, axis)
	log("Sliced %d entries down to %d data points" % (fulldata.size, sliced.size))

	np.savetxt(
		fname=outputfh,
		X=sliced,
		fmt=args.nformat,
		header=','.join(['"%s"' % s for s in sliced.dtype.names]),
		comments="# %s\n" % description.replace('slices', 'sliced')
	)
