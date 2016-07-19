# presence of this file defines python moudule `exaiohelpers`
from __future__ import print_function # Py3 compilance
from functools import wraps
import sys

# convenient python2 methods

verbose = False

def log(text, newline=True, force=True):
	if not force and not verbose: return
	print(text, file=sys.stderr, end="\n" if newline else '')
	sys.stderr.flush()

def vectorize_concatenate(func):
	"""
	Decorator to vectorize file input readers like given in argio or vtkreader.
	"""
	import numpy as np
	def func_wrapper(*args, **kwargs):
		args_without_fname = list(args)
		fnames = args_without_fname.pop(0)
		if not isinstance(fnames ,basestring):
			# progress bar: works but not nice output
			#genlogwrap = lambda i: lambda text,*a,**kw: log("[%02d/%02d] "% (i,len(fname))+text,*a,**kw)
			
			# optimization to avoid np.vstack for large input files (>3GB):
			if len(fnames) == 1:
				return func(fnames[0], *args_without_fname, **kwargs)

			outputs = [func(f, *args_without_fname, **kwargs) for i,f in enumerate(fnames)]
			return np.vstack(tuple(outputs))
		else:
			return func(*args, **kwargs)
	return func_wrapper