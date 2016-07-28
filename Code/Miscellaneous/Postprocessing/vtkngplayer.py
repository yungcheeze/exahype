#!/usr/bin/env python
#
# This script is part of ExaHyPE. Probably.
# This script is Python 2.
#
# A minimalistic interactive movie player for 1D or 2D data
#
# -- Public Domain, Jul, 22. 2016, Sven K.

from __future__ import print_function
from exaiohelpers import argio

import numpy as np
import matplotlib.pyplot as plt

# Python built in batteries
import argparse, sys # Python 2.7
import gzip

programname = sys.argv[0]
parser = argparse.ArgumentParser(description='ExaHyPE simulation data player',
	epilog="\n".join([
		"A minimalistic interactive movie player for 1D or 2D data.",
		"currently only displays x-y layer or x axis."
	]), formatter_class=argparse.RawDescriptionHelpFormatter)


argio.add_io_group(parser)
args = parser.parse_args()
log = argio.logger_for(args)
sol = argio.get_input(args)

# TODO: Should get 1d or 2d parameters from command line
dim = 2

# TODO: should be removed at some time for noninteractive sessions
plt.ion();

def find_nearest(a, a0):
	return a.flat[abs(a-a0).argmin()]
unzip = lambda z: zip(*z) # inverse zip

uni = lambda col: unique(sol[col])

times = uni('time')
x, y, z = map(uni, "xyz")
spacedim = 3 if len(z) > 1 else 2

log("vtkplayer: Having loaded %d timesteps in the range (%f,%f) in a %d-dimensional problem" % (len(times), min(times), max(times), spacedim))

# Quantities: name of conserved variable in ExaHypE -[mappedTo]-> human readable name
# TODO: Should be get from configfile or so
quantities = ["Q0", "Q1", "Q4"]
quantities_readable = ["d", "sx", "e"]
yslice = y[int(len(y)/2)]
zslice = z[0]

def slice2d(time, quantity):
	data2d = sol[ (sol['time'] == time) & (sol['z'] == zslice) ]
	X, Y = meshgrid(x, y)

	# stopped at this point and reported bug #50.
	
	#sx, ix = unique(data2d['x'], return_index=True)
	#sy, iy = unique(data2d['y'], return_index=True)
	#d = data2d[quantity][
	#sy = data1d[quantity][ix]
	return (sx,sy)

def slice1d(time, quantity):
	data1d = sol[ (sol['time'] == time) & (sol['y'] == yslice) & (sol['z'] == zslice) ]
	sx, ix = unique(data1d['x'], return_index=True)
	sy = data1d[quantity][ix]
	return (sx,sy)

def createPlot(ax, quantity):
	sx, sy = slice1d(times[0], quantity=quantity)
	line, = ax.plot(sx, sy, "o-")
	title = ax.set_title("title")
	ax.set_xlim(0,1) # ExaHyPE coordinate system range
	return (line,title)

def draw_update(time_slider=0):
	time = find_nearest(times, time_slider)
	print "Plot for time="+str(time)
	
	for i, qu in enumerate(quantities):
		tt = "Slice of %s at t=%f, y=%.1f, z=%.1f" % (quantities_readable[i], time, yslice, zslice )
		axes[i].set_title(tt)
		lines[i].set_data(slice1d(time, quantity=qu))
		show()

gs = GridSpec(len(quantities),1)
axes = map(subplot, gs)
lines, titles = unzip([ createPlot(a,q) for a,q in zip(axes, quantities)])
subplots_adjust(left=0.10, top=0.85)
axcolor = 'lightgoldenrodyellow'
ax1 = plt.axes([0.10, 0.90, 0.8, 0.03], axisbg=axcolor)

timeslider = Slider(ax1, "Time", min(times), max(times), min(times))
timeslider.on_changed(draw_update)

draw()


### TODO for the 2D movies: Create multiple panels which show all quantities at the same
### time, next to each other. So basically the same as for 2D just flipped due to 16:9.
