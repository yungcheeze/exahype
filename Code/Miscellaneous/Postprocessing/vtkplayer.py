#!/usr/bin/python
#
# A very simple Python2 tool to introspect 1D slices from 2D (3D)
# ExaHyPE simulation data. Used so far for the shocktube.
#
# -- Sven K, 2016

import sys
from pylab import *
ion();

def find_nearest(a, a0):
	return a.flat[abs(a-a0).argmin()]
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

print "vtkplayer: Having loaded %d timesteps in the range (%f,%f) in a %d-dimensional problem" % (len(times), min(times), max(times), dim)

# Quantities: name of conserved variable in ExaHypE -[mappedTo]-> human readable name
quantities = ["Q0", "Q1", "Q3"]
quantities_readable = ["d", "sx", "e"]
yslice = y[int(len(y)/2)]
zslice = z[0]

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

