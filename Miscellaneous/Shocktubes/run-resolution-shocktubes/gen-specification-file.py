#!/usr/bin/env python

from __future__ import print_function
import re
from numpy import arange, transpose

# Peano divides each refinement level by 3
peanothird = 1./3.

# number of levels wanted
level = arange(1,7)

dx = peanothird ** level
dxeps = 1.2 * dx
num_points = 1./dx

print("Requested levels=", level)
print("Giving dx=", dx)
print("Will print out slightly larger dxeps=", dxeps)
print("And number of grid points=", num_points)

template_file = "SRHD.exahype.tpl"
target_filename = "SRHD-l%d.exahype"
template = open(template_file, 'r').read()

def create_tpl_for(cur_level):
	level_index = list(level).index(cur_level)
	def repl(m):
		key = m.group(1).lower()
		try:
			replacement_array = globals()[key]
			value = replacement_array[level_index]
			return str(value)
		except KeyError:
			print("Template variable @%s@ does not exist" % key)

	configfile = re.sub('@([a-z_A-Z]+)@', repl, template)
	filename = target_filename % cur_level
	with open(filename, 'w') as outfile:
		outfile.write("// GENERATED CONFIGURATION FILE\n")
		outfile.write(configfile)
	print("Finished generating ", filename)

#level_indices = transpose(list(enumerate(level)))[0] # if there would be a nicer way...
map(create_tpl_for, level)

print("Done")
