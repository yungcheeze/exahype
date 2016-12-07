#!/usr/bin/env python
#
# This is a python2 module with helper functions for doing convergence studies
# with ExaHyPE. Should also work with python3.
#
# SK, 2016

import os, subprocess, sys
from os import path, stat
from functools import partial

# provide StringIO consistently
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO

shell = lambda cmd: subprocess.check_output(cmd, shell=True).strip()
getenv = lambda key, default=None: os.environ[key] if key in os.environ else default
pidlist = lambda processes: " ".join([str(proc.pid) for proc in processes])

def runBinary(binary, envUpdate):
	"""
	Runs the command "binary", this may be only one program without parameters,
	with the environment which is composed of the environmental variables and
	the envUpdate dictionary.
	Example usage:
	> runBinary("ExaHyPE-FooBar", { 'EXAHYPE_SKIP_TEST': 'True' })

	You can replace subprocess.Popen with other alternatives, but make sure
	they return some kind of process handle which is subsequently used
	"""
	env = os.environ.copy()
	env.update(envUpdate)
	
	if not os.path.exists(binary):
		raise IOError("Failure: '%s' does not exist" % binary)

	return subprocess.Popen([binary], env=env)

def read_simulation_params(envfile):
	"""
	Small function to parse files in this format:
	EXAMESHSIZE=0.0037037037037
	EXAPORDER=2.0
	EXAHYPE_INITIALDATA=MovingGauss2D
	...
	Returns dictionary like
	{
	'EXAMESHSIZE': 0.0037037037037,
	'EXAPORDER': 2.0,
	'EXAHYPE_INITIALDATA': 'MovingGauss2D',
	...
	}
	"""
	from re import match, sub
	params = {}
	with open(envfile, 'r') as fh:
		for line in fh:
			res = match(r'^([a-zA-Z0-9_]+)=(.*)$', line)
			if res:
				k, v = res.group(1), res.group(2)
				# try to cast the value as a float,
				# as this is a common data type in our buisness
				try:
					params[k] = float(v)
				except:
					params[k] = v
	return params

is_empty_file = lambda f: stat(f).st_size == 0
simfile = lambda fname: lambda simdir: path.join(simdir, fname)

def trygenfromtxt(f, *args, **kwargs):
	"fault tolerant genfromtxt"
	try:
		return np.genfromtxt(f, *args, **kwargs)
	except:
		return None

def is_headless():
	return not 'DISPLAY' in os.environ

def gensvg(figure):
	"""
	Save a matplotlib figure to an SVG (scalable vector graphics) file
	and returns the file contents as string.
	"""
	#import matplotlib
	#matplotlib.use('Agg') # if headless
	#import matplotlib.pyplot as plt
	import StringIO

	imgdata = StringIO.StringIO()
	figure.savefig(imgdata, format='svg')
	imgdata.seek(0)
	return imgdata.buf

def shortenPathsInTable(dftable, possiblecolumns):
	for col in possiblecolumns:
		try: dftable[col] = dftable[col].map(path.basename)
		except e: pass
	
# another pandas trick to remove constant columns
# and one to retrieve the stripped values
def stripConstantColumns(dftable):
	"""
	Detects columns with constant values in a Pandas Dataframe and returns
	a tuple (A,B) where A is a dataframe *without* the constant columns and
	B is a dataframe with all the constant column values transposed.
	"""
	import pandas as pd
	reducedtable = dftable.loc[:, (dftable != dftable.ix[0]).any()]
	removedcolumns = set(dftable.columns) - set(reducedtable.columns)
	removedvalues = {c: dftable.ix[0][c] for c in removedcolumns}
	removedtable = pd.DataFrame({ 'Column names': removedvalues.keys(), 'Constant values': removedvalues.values() })
	assert len(dftable.columns) - len(reducedtable.columns) == len(removedtable)
	return (reducedtable, removedtable)

def keepColumnsIf(pdtbl, testColumNameFunction):
	"""
	Keep columns in pandas table where testColumNameFunction(c) gives
	true on the column name. This is a more generic version of pdtbl.filter()
	which only supports wildcards and regexp.
	"""
	cols = [col for col in pdtbl.columns if testColumNameFunction(col)]
	return pdtbl[cols]

def MapColumnNames(pdttbl, ColumnMapFunction, inplace=True):
	"""
	Allows mapping names of columns in pandas table.
	Can be used inplace or not. In any case returns new pdttbl.
	"""
	cols = list(pdttbl.columns)
	newcols = map(ColumnMapFunction, cols)
	return pdttbl.rename(columns=dict(zip(cols, newcols)), inplace=inplace)

def RemoveStringInColumns(pdttbl, removalstring, inplace=True):
	"Shorthand when you just want to remove stuff"
	return MapColumnNames(pdttbl, lambda k: k.replace(removalstring, ''), inplace)


def parse_time(regex, time_str):
	# cf http://stackoverflow.com/a/4628148
	from datetime import timedelta
	import re
	parts = re.compile(regex).match(time_str)
	if not parts:
		return
	parts = parts.groupdict()
	time_params = {}
	for (name, param) in parts.iteritems():
		if param:
			time_params[name] = int(param)
	return timedelta(**time_params)

# gives a function to parse time on
time_parser = lambda regex: partial(parse_time, regex)
	
def executeTemplate(inputfilename, tmplvars, outputfilename):
	"""
	Interpretes a file as python template string, inserts the tmplvars and
	writes output to outputfile. Usage example: Let content of inputfile be
		<p>This is a convergence report generated at <strong>${DATE}</strong> on ${HOST} by ${WHOAMI}.
	And your tmplvars be a dictionary
		{ 'DATE': 'today', 'HOST': 'me', 'WHOAMI': 'notyou' }
	Then the resulting string which goes to outputfile is obvious, isn't it?
	"""
	import string
	tmpl=open(inputfilename, 'r').read().strip()
	html = string.Template(tmpl).substitute(tmplvars)
	with open(outputfilename, 'w') as out:
		out.write(html)

class Template:
	"""
	A classy interface to executeTemplate for delayed execution. Usage:
	> tpl = Template('./my-template.html', './my-output.html')
	> tpl.set('something_global', 'foo') # can be done before
	> tpl.execute({ 'template': 'vars', 'go': 'here' }) # or at runtime
	"""
	def __init__(self, inputfile, outputfile):
		self.inputfile = inputfile
		self.outputfile = outputfile
		self.prepared_args = {}
	def set(self, key, value):
		self.prepared_args[key] = value
	def execute(self, tmplvars=dict(), verbose=False):
		self.prepared_args.update(tmplvars)
		executeTemplate(self.inputfile, self.prepared_args, self.outputfile)
		if verbose:
			print "Wrote report to %s." % self.outputfile
	def addStatistics(self):
		from datetime import datetime
		from socket import gethostname
		from getpass import getuser

		self.set('DATE', datetime.now().strftime("%c"))
		self.set('HOST', gethostname())
		self.set('WHOAMI', getuser())
		return self # chainable
