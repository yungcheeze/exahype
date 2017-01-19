#!/usr/bin/env python
#
# This is a standalone script to produce a two-panel figure where errorplots
# and timestep plots are just above each other.
#
# Input sources are the ASCII reduction integrals (error norms) as well as
# the exahype log file containing the dt_min, etc. information
#
#
# To be used from ipython, in the moment.


import matplotlib.pyplot as plt
import extract_dtmin
from os import path
import pandas as pd
from glob import glob
import sys

plt.ion()
possiblecolors = "r g k b m y r c m".split(" ")

def timesteps(logfile):
	"Reads logfile (-> csv) -> dataframe"
	print "Extracting timesteps from %s..." % logfile
	m = extract_dtmin.memorywriter()
	extract_dtmin.timesteps(logfile, m.writer)
	return pd.DataFrame(m.storage, dtype=float)

def errorplots(simulationsdirs):
	#errorPlot = plt.figure(figsize=(18,8))
	xcolumn = 'time'
	ycolumn = 'max' # l2norm, l1norm
	outputdir = 'output'
	quantity = 'error-rho.asc'

	num_markers = 10
	for i,simulationdir in enumerate(simulationsdirs):
		quantityfile = path.join(simulationdir, outputdir, quantity)
		data = pd.read_csv(quantityfile, delim_whitespace=True)
		plt.plot(data[xcolumn], data[ycolumn], label=simulationdir, color=possiblecolors[i])

	l = plt.legend(loc='center left')#, bbox_to_anchor=(1,0.5))
	l.draggable()
	#plt.subplots_adjust(right=0.8) # was space for legend, no more.
	plt.title("Error evolution: %s of quantity %s" % (ycolumn, quantity))
	plt.xlabel("Simulation time")
	plt.ylabel("Error")
	ax = plt.gca()
	ax.set_yscale('log')
	plt.ylim(1e-10, 1e0)

def dtplots(simulationdirs):
	for i, simulationdir in enumerate(simulationdirs):
		logfile = glob(path.join(simulationdir, '*.log'))[0]
		data = timesteps(logfile)
		plt.plot(data['t_min'], data['dt_min'], label=simulationdir, color=possiblecolors[i])
	plt.xlabel("Simulation time")
	plt.ylabel("Time step")

def combinedplots(simulationdirs):
	plt.subplot(211)
	errorplots(simulationdirs)
	plt.subplot(212)
	dtplots(simulationdirs)

if __name__ == "__main__":
	combinedplots(sys.argv)

