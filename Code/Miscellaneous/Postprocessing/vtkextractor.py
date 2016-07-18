#!/usr/bin/python
#
# This script is part of ExaHyPE. Probably.
# This script is Python 2.
#
# This is a small helper script to convert VTK ExaHype files to plain ASCII. These
# are typically a set of files "solution-0.vtk, solution-1.vtk, ... solution-50.vtk, ..."
# with each 200MB in size, containing vtkUnstructuredGrids. 
# This script collects them and produces a single output file which may be ASCII
# for convenience.
#
# This script is not performant at all, it collects all input in numpy arrays before
# writing any output. In a typical scenario with 40 output files, 3GB in size, the
# script needs around 8 mins and 2GB system memory before writing a 200MB compressed
# output file.
#
# -- Public Domain, May, 10. 2016, Sven K.

from __future__ import print_function
from exaiohelpers import argio
import numpy as np

# Python built in batteries
import argparse, sys # Python 2.7
import gzip

programname = sys.argv[0]
parser = argparse.ArgumentParser(description='ExaHyPE simulation data converter, ie. VTK to plain ASCII converter',
	epilog="\n".join(["Example invocations:",
		"   %s -o solution-0.txt.gz -c solution-0.vtk" % programname,
		"   %s solution-*.vtk | gzip > foobar.csv.gz" % programname,
		"if you have files solution-1.vtk ... solution-15.vtk ... than you probably",
		"want natural sorting, eg. with",
		"   %s -o data.csv.gz -c $(ls *vtk | sort -V)" % programname,
		"%s will not never change the order of the input data." % programname
	]),
	formatter_class=argparse.RawDescriptionHelpFormatter)


argio.add_io_group(parser)
args = parser.parse_args()
log = argio.logger_for(args)

data = argio.get_input(args)
log("Have read a %s-shaped numpy array", data.shape)
argio.write_output(data, args)
log("Finished")

