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
from vtkreader import vtkreader, log
import numpy as np

# Python built in batteries
import argparse, sys # Python 2.7
import gzip

output_formats = {}
def register_format(command_name):
	def decorator(func):
		output_formats[command_name] = func	
		return func
	return decorator

@register_format("np")
def output_npy(npdata, outputfname, args=None):
	np.save(outputfname, npdata)

@register_format("csv")
def output_ascii(npdata, outputfname, args):
	log("This may take a while. Have a coffee. Or watch your output growing.")
	# TODO: I inserted a bug here. Doesn't work any more with ASCII here. Dunno why.
	np.savetxt(
		fname=outputfname,
		X=npdata,
		fmt=args.nformat,
		header=','.join(['"%s"' % s for s in npdata.dtype.names]),
		comments=''
	)


programname = sys.argv[0]
parser = argparse.ArgumentParser(description='ExaHyPE VTK to plain ASCII converter',
	epilog="\n".join(["Example invocations:",
		"   %s -o solution-0.txt.gz -c solution-0.vtk" % programname,
		"   %s solution-*.vtk | gzip > foobar.csv.gz" % programname,
		"if you have files solution-1.vtk ... solution-15.vtk ... than you probably",
		"want natural sorting, eg. with",
		"   %s -o data.csv.gz -c $(ls *vtk | sort -V)" % programname,
		"%s will not never change the order of the input data." % programname
	]),
	formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('fnames', metavar='solution-0.vtk', nargs='+',
                   help='The VTK file(s) to convert. Multiple input files are simply attached in output.')
parser.add_argument('-c', '--compress', dest='compress', action='store_true', default=False,
                   help='Compress the output, creates a gzip file')
parser.add_argument('-q', '--quiet', dest='quiet', action='store_true', default=False,
                   help='Quiet, no informal output on stderr')
parser.add_argument('-o', '--output', dest='output', default=False,
                   help='Output filename (default: stdout)')
parser.add_argument('-f', '--format', dest='format', choices=output_formats.keys(), type=str, default='csv',
		   help='Output file formats. CSV (ASCII) takes long, binary formats are rather quick.')
parser.add_argument('-n', '--numberformat', dest='nformat', default='%.5e',
		   help='Number format string, cf. numpys savetxt() documentation. Applies only for CSV.')


args = parser.parse_args()
use_stdout = not args.output

if args.quiet:
	def log(text, newline=True, force=True):
		pass

opener = gzip.open if args.compress else open
if args.compress and use_stdout:
	log("Gzip output for stdout currently not supported. Just use pipes: %s ... | gzip" % programname, force=True)
	sys.exit(1)
outputfh = sys.stdout if not args.output else opener(args.output, 'w')

npdata = vtkreader(args.fnames, log)
log("Printing output to %s" % str(outputfh))

try:
	output_formats[args.format](npdata, outputfh, args)
except IOError:
	if use_stdout:
		# eg. when using this command | head
		# stdout is closed, no point in continuing
		# Attempt to close them explicitly to prevent cleanup problems:
		try:
			sys.stdout.close()
		except IOError:
			pass
		try:
			sys.stderr.close()
		except IOError:
			pass
log("Finished")

