#!/usr/bin/python
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
# -- Public Domain, 2016, Sven K.

from __future__ import print_function # Py3 compilance

from vtk import vtkUnstructuredGridReader
from vtk.util import numpy_support
import numpy as np

# Python built in batteries
import argparse, sys # Python 2.7
from collections import namedtuple # Python 2.4
from functools import wraps
import gzip

verbose = False

def log(text, newline=True, force=True):
	if not force and not verbose: return
	print(text, file=sys.stderr, end="\n" if newline else '')
	sys.stderr.flush()

def vtkreader(fname, log=log):
	"""
	Transforms a single VTK files to a recorded numpy array.
	Currently only filenames given as string are supported, no python streams.
	For a 200MB input file, this function needs around 1-2sec.
	"""

	# "vectorize" this function for list of fnames
	if not isinstance(fname ,basestring):
		# works but not nice output
		#genlogwrap = lambda i: lambda text,*a,**kw: log("[%02d/%02d] "% (i,len(fname))+text,*a,**kw)
		outputs = [vtkreader(f, log) for i,f in enumerate(fname)]
		return np.vstack(tuple(outputs))
	
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

	ret = datasets.transpose()
	#numpy.float32
	structured_array_dtype = np.dtype([ (k, ret.dtype) for k in column_names ])
	ret.dtype = structured_array_dtype
	return ret


