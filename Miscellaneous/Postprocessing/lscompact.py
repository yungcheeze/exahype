#!/usr/bin/python3
# a compact directory listing
# (unix) command line tool
# Pure Python3 (no dependencies)
# (GPL 2016) SvenK

from sys import exit
from os import listdir, path
from re import sub, search, split
from functools import partial
from itertools import groupby

# when `rgxp` is a regular expression, ensure it matches on full strings
fullRegex = lambda rgxp: "^" + rgxp + "$"

removeDuplicates = lambda lst: list(set(lst))
applyValues = lambda fn, dct: { k: fn(v) for k,v in dct.items() }
mapValues = lambda fn, dct: applyValues(lambda v: list(map(fn, v)), dct)
nest = lambda f, g: lambda *a, **kw: f(g(*a,**kw))
flatten = lambda lst: [item for sublist in lst for item in sublist]
complement = lambda first, second: [x for x in first if x not in second]
at = lambda L, I: [ L[i] for i in I ]

def unique(lst, return_index=True):
	"Unique sublist of `lst` returning index array as in numpy.unique"
	I, U = [], []
	for i, e in enumerate(lst):
		if not e in U:
			U.append(e)
			I.append(i)
	return (U,I) if return_index else U

def naturalSort(lst):
	"human/natural sorting of a list (ie. 1, 10, 100 like 001, 010, 100)"
	# source: http://stackoverflow.com/a/4836734
	convert = lambda text: int(text) if text.isdigit() else text.lower() 
	alphanum_key = lambda key: [ convert(c) for c in split('([0-9]+)', key) ] 
	return sorted(lst, key = alphanum_key)

def ranges(lst):
	"Identify groups of continous numbers in a list"
	# source: http://stackoverflow.com/a/2154409
	pos = (j - i for i, j in enumerate(lst))
	t = 0
	for i, els in groupby(pos):
		l = len(list(els))
		el = lst[t]
		t += l
		yield (el, el+l)


def lscompact(allfiles, placeholder="<>", numdetector=r'\d+', expandable=False, recursive=False):
	"""
	Compact listing of directories.
	"""
	expansion_pattern = "[%d-%d]"
	if expandable: expansion_pattern = "{%d,%d}"

	# seperate interesting filenames with numbers from those without
	numberfiles = list(filter(partial(search, numdetector), allfiles))
	nonnumericfiles = complement(allfiles, numberfiles)
	# patterns makes all solution3.vtk to solution<>.vtk
	patterns = list(map(partial(sub, numdetector, placeholder, count=1), numberfiles))
	# regexes contains lists which allow filtering similiar strings
	regexes = list(map(nest(fullRegex, partial(sub, numdetector, numdetector)), numberfiles))
	# we do the filtering by removing similar regexes
	uniqueRegexes, uniqueIdx = unique(regexes)
	patterns = at(patterns, uniqueIdx)
	assert len(patterns) == len(uniqueRegexes)
	# find the matching files for each matching regexp
	matching = { pattern: naturalSort(filter(partial(search, regexp), numberfiles)) for regexp, pattern in zip(uniqueRegexes, patterns) }
	# determine the numbers in these filenames
	numbers = mapValues(nest(int, partial(sub, r'.*?(\d+).*', r'\1')), matching)
	# find groups of numbers in the filenames
	groups = applyValues(nest(list, ranges), numbers)
	# replace the group(start, end+1) objects with filename representatives
	speakingRange = lambda template: lambda rng: template.replace(placeholder, expansion_pattern % (rng[0], rng[1]-1) if rng[1]-rng[0] > 1 else str(rng[0]))
	speakingGroups = flatten([ list(map(speakingRange(template), ranges)) for template, ranges in groups.items() ])
	# merge the reduced numeric files and the nonnumeric files
	reducedfiles = sorted(speakingGroups + nonnumericfiles)

	# currently: poor man's print
	for f in reducedfiles:
		print(f)
		if(path.isdir(f)):
			# this is not really supported yet:
			lscompact(listdir(f), placeholder=placeholder, numdetector=numdetector, expandable=expandable, recursive=recursive)

	# ideas for better output:

	# not adaptive, contains ['python lists markup']
	# from pprint import pprint
	# pprint(reducedfiles, compact=True)

	# adaptive multicolumn like /bin/ls could be archived with
	# https://gist.github.com/jtriley/1108174
	# and http://stackoverflow.com/a/25027568

def cli():
	import argparse
	parser = argparse.ArgumentParser(description="Compact directory listing")
	parser.add_argument('-e', '--expandable', action='store_true', help='Write names shell-expandable like image{01,10}.png')
	parser.add_argument('-r', '--recursive', action='store_true', help='Go recursively into subdirectories (like find)')
	parser.add_argument('files', nargs='*', help='List of files or directories, as for /bin/ls')
	args = parser.parse_args()
	# if no arguments given, list the current directory.
	files = args.files if len(args.files) else listdir(".")
	lscompact(files, expandable=args.expandable, recursive=args.recursive)

if __name__ == "__main__":
	cli()


