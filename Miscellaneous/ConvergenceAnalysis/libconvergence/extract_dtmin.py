#!/usr/bin/python

"""
This script collects information about a timestep, as in

2017-01-19 14:15:55  0.770825     info         step 14  t_min          =0.158573
2017-01-19 14:15:55  0.770835     info          dt_min         =0.0113432
2017-01-19 14:15:55  0.77086      info          memoryUsage    =29 MB

together to an ASCII file as in 

step	t_min		dt_min
14	0.158573	0.0113432


To invoke, call it on a logfile of an ExaHyPE run,
or use as library

"""

import sys, re, csv

keys = ["step", "t_min", "dt_min"]

class csvwriter:
	def __init__(self, outfile=sys.stdout):
		self.csvwriter = csv.writer(sys.stdout, delimiter="\t")
		self.csvwriter.writerow(keys) # writes csv header

	def writer(self, outrow):
		self.csvwriter.writerow([ outrow[k] for k in keys ])

class memorywriter:
	def __init__(self):
		self.storage = []
	def writer(self, outrow):
		self.storage.append(outrow)

def timesteps(exahype_logfile, writer):
	outrow = {}
	fh = open(exahype_logfile, 'r')

	for line in fh:
		if "info" in line:
			stepline = re.match(r".*step\s+(?P<step>\d+).*", line)
			if stepline:
				if len(outrow):
					writer(outrow)
				outrow = {}
				outrow['step'] = stepline.group('step')
				t_min = re.match(r".*t_min\s+=(?P<t_min>[\d\.]+)\s*", line)
				if t_min:
					outrow['t_min'] = t_min.group('t_min')
			else:
				dt_min = re.match(r".*dt_min\s+=(?P<dt_min>[\d\.]+)\s*", line)
				if dt_min:
					outrow['dt_min'] = dt_min.group('dt_min')

# convenient functions:

def timesteps2csv(logfile, outfile=sys.stdout):
	"Reads logfile -> csv"
	c = csvwriter(outfile)
	timesteps(logfile, c.writer)

def timesteps2dataframe(logfile):
	"Reads logfile -> dataframe"
	import pandas as pd
	m = memorywriter()
	timesteps(logfile, m.writer)
	return pd.DataFrame(m.storage, dtype=float)

if __name__ == "__main__":
	if not len(sys.argv) == 2:
		print "Usage: %s <exahype-log-file>" % sys.argv[0]
		sys.exit(-1)
	fname = sys.argv[1]
	timesteps2csv(fname)

