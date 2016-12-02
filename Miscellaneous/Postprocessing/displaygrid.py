#!/usr/bin/python env
#
# Display the gridpoints. Only useful for very small grid resolvement levels,
# eg. L=3, with 3**L = 27 grid points. Beyond that the output of this gets crazy.
#

from pylab import *
ion()

sol = load("solutions.np")

uni = lambda col: unique(sol[col])
middle = lambda a: a[int(len(a)/2)]

times = uni('time')
x, y, z = map(uni, "xyz")
time = middle(times)
zslice = middle(z)

data2d = sol[ (sol['time'] == time) & (sol['z'] == zslice) ]
X, Y = data2d['x'], data2d['y']

# X and Y now holds all the (matching) X and Y coordinates
# of the individual data points in the solutions. At each
# data point, ALL quantities are available.

all_points, = plot(X, Y, "x-", color="lightgrey")

# previous points (initial values do not matter)
hist, = plot([0],[0], "bo-")
prev, = plot([0],[0], "rD-")

title("Peano writes two times too much data")
xlabel("x axis")
ylabel("y axis")

def update_plot(up):
	if up > 1:
		# show the last points
		hist.set_data(X[:up-1], Y[:up-1])
		prev.set_data(X[up-1], Y[up-1])

	#main.set_data(X[:up], Y[:up])
	draw()
	show()

# make an animation
do_animation = True
output_filename = "animation.mpg"
if do_animation:
	import matplotlib.animation as manimation
	FFMpegWriter = manimation.writers['ffmpeg']
	writer = FFMpegWriter(fps=10)
	with writer.saving(gcf(), output_filename, 100):
		for i,_ in enumerate(X[:200]):
			print "Movie: %d from %d slides done." % (i, len(X))
			update_plot(i)
			writer.grab_frame()
