#!/usr/bin/python
#
# A very simple Python2 tool give out an ASCII table about the values of
# the simulation on a 1D slice for a given time. Very tailored for the shocktube.
#
# As the first argument, give the path to numpy file.
#
# As the second argument, you have give the time index you are interested in. Example
# values are "0" for the very first time step or "-1" for the last time step.
# Without an argument, all time steps will be given.
#
# -- Sven K, 2016

import sys
from pylab import *
ion();

unzip = lambda z: zip(*z) # inverse zip

print sys.argv

la = len(sys.argv)

if la not in [2,3]:
	print "Usage: %s <path to numpy file> <timeindex>" % sys.argv[0]
	print "Read my code for further documentation"
	sys.exit(-1)

times_info_mode = (la == 2)
timeindex = int(sys.argv[2]) if not times_info_mode else 0

fname = sys.argv[1]
sol = load(fname)
uni = lambda col: unique(sol[col])

times = uni('time')
x, y, z = map(uni, "xyz")
dim = 3 if len(z) > 1 else 2

print >>sys.stderr, "vtkfinal1d: Having loaded %d timesteps in the range (%f,%f) in a %d-dimensional problem" % (len(times), min(times), max(times), dim)

# Quantities: name of conserved variable in ExaHypE -[mappedTo]-> human readable name
quantities = ["Q0", "Q1", "Q4"]
quantities_readable = ["d", "sx", "e"]
yslice = y[int(len(y)/2)]
zslice = z[0]

if times_info_mode:
	print >>sys.stderr, "Printing information about times"
	print "it t"
	for it, t in enumerate(times):
		print "%d %f" % (it, t)
	sys.exit(0)

print >>sys.stderr, "Printing for time index = %d" % timeindex

try:
	time = times[timeindex] # last snapshot
except:
	print "Please give a correct time. Run without arguments to see valid time indices"
	sys.exit(-1)

print >>sys.stderr, "Printing for time = %f" % time

def slice1d(quantity):
	data1d = sol[ (sol['time'] == time) & (sol['y'] == yslice) & (sol['z'] == zslice) ]
	sx, ix = unique(data1d['x'], return_index=True)
	return data1d[quantity][ix]

columns = map(slice1d, quantities)
columns.insert(0, x)
colnames = "x d sx e"

savetxt(sys.stdout, transpose(columns), fmt='%.8e', header=colnames, comments='')

