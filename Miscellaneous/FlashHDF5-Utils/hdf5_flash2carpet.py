#!/usr/bin/python
# Python 2. Needs HDF5 for Python (h5py).
#
# This script translates a FlashHDF5 file (as created by ExaHyPE) to
# a CarpetHDF5 file. This file can then read by the usual readers which
# we have, for instance in any Visit installation.
#
# Current limitations:
# * Not tested for multiple MPI files. Will probably just work via
#   translating each file on its own.
#
# -- SvenK, 2017-06-12
#

import sys
import numpy, h5py

if (not len(sys.argv[1])):
	print "Usage: flash2carpet.py <InputFlashHDF5File.h5> <OutputCarpetHDF5File.h5>"

flash_filename = sys.argv[1]
carpet_filename = sys.argv[2]

print "Reading from ", flash_filename
print "Writing to ", carpet_filename

flash = h5py.File(flash_filename, "r")
carpet = h5py.File(carpet_filename, "w")

# Write carpet global group
ranks = 1
g = carpet.create_group("Parameters and Global Attributes")
g.attrs.create("nioprocs", ranks, dtype="int32")


for groupname, group in flash.iteritems(): # ie. go over iteration/time
	if not "fields" in group:
		print "Skipping group", groupname
		continue

	fields = group["fields"]
	components = fields.shape[0]
	patches = fields.shape[1:-1]
	dim = len(patches) # 1, 2 or 3 dimensional
	patchslice = [slice(None)] * dim
	unknowns = fields.shape[-1] # number of written unknowns
	time = group.attrs["time"] # timestamp
	it = group.attrs["iteration"]

	for comp in range(components):
		origin = group["origin"][comp]
		dx = group["delta"][comp]
		assert origin.shape == (dim,), "Dimensionality mismatch"
		for unknown in range(unknowns):
			# currently missing: a correct dataset name
			name = "ExaHyPE::flash%d" % unknown
			data = fields[ tuple([comp] + patchslice + [unknown]) ]
			assert data.shape == patches, "Somehow data are misformed"
			dset = carpet.create_dataset("%s it=%d tl=0 m=0 rl=0 c=%d" % (name,it,comp),
				shape=data.shape, dtype=data.dtype, data=data)

			dset.attrs["origin"] = origin
			dset.attrs["iorigin"] = numpy.array([0]*dim,dtype='int32')
			dset.attrs.create("level", 0, dtype='int32')
			dset.attrs.create("timestep", it, dtype='int32')
			dset.attrs.create("time", time, dtype='float')
			dset.attrs["delta"] = numpy.array(dx, dtype='float')
			dset.attrs["name"] = numpy.string_(name + "\0")

