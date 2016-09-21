#!/usr/bin/env python

from numpy import genfromtxt
# batteries:
from glob import glob
from os import path

# which quantity shall we look at, in the moment?
quantityfile="output/error-rho.asc"

simbase = "simulations" # simulation directory

# genfromtxt with header detection is very sensible to the first line's format.
# this will not work:
#    # 1:plotindex ,2:time ,3:l1norm ,4:l2norm ,5:max ,6:min ,7:avg ,
# in constrast, the line has to look like
#    1:plotindex 2:time 3:l1norm 4:l2norm 5:max 6:min 7:avg
# or even better
#    plotindex time l1norm l2norm max min avg

files = glob(path.join(simbase, '*', quantityfile))

[genfromtxt(f, names=True) for f in files]
