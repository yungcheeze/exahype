#!/usr/bin/python
#
# A very simple Python2 tool give out an ASCII table about the final
# values of the simulation on a 1D slice. Very tailored for the shocktube.
#
# -- Sven K, 2016

import sys
from pylab import *
ion();

unzip = lambda z: zip(*z) # inverse zip

if len(sys.argv) != 2:
	print "Usage: %0 <path to numpy file>" % sys.argv[0]
	sys.exit(-1)

fname = sys.argv[1]
sol = load(fname)
uni = lambda col: unique(sol[col])

times = uni('time')
x, y, z = map(uni, "xyz")
dim = 3 if len(z) > 1 else 2

print "vtkfinal1d: Having loaded %d timesteps in the range (%f,%f) in a %d-dimensional problem" % (len(times), min(times), max(times), dim)

# Quantities: name of conserved variable in ExaHypE -[mappedTo]-> human readable name
quantities = ["Q0", "Q1", "Q3"]
quantities_readable = ["d", "sx", "e"]
yslice = y[int(len(y)/2)]
zslice = z[0]

def slice1d(quantity):
	time = sol['time'][-1] # last snapshot
	data1d = sol[ (sol['time'] == time) & (sol['y'] == yslice) & (sol['z'] == zslice) ]
	sx, ix = unique(data1d['x'], return_index=True)
	return data1d[quantity][ix]

columns = map(slice1d, quantities)
columns.insert(0, x)
colnames = "x d sx e"

savetxt(sys.stdout, transpose(columns), fmt='%.8e', header=colnames)

