#!/usr/bin/python
"""
.. module:: main
  :platform: Unix, Windows, Mac
  :synopsis: Generates reference basis functions and their first and second derivatives up to a specific N.

:synopsis: Generates reference basis functions and their first and second derivatives up to a specific N.
"""
from lagrangeinterp import *
import numpy as np
import sympy
from sympy.printing import print_ccode
import re

max_order = 9
# Print the 1d reference basis functions#
print "// Basis functions in interval (0,1)."
for N in range(0,max_order+1):
    for m in range(0,N+1):
        print("/** Lagrange basis polynomial of order %d which passes through the %d-th support point in interval (0,1). */" % (N,m))
        print("double basisFunction_%d_%d(const double s);" % (N,m))
    print ""

# Print the first derivative of the 1d reference basis functions
print "\n\n// First derivatives of basis functions in interval (0,1)."
for N in range(0,max_order+1):
    for m in range(0,N+1):
        print("/** First derivative of Lagrange basis polynomial of order %d which passes through the %d-th support point in interval (0,1). */" % (N,m))
        print("double basisFunctionFirstDerivative_%d_%d(const double s);"  % (N,m))
    print ""

## Print the second derivative of the 1d reference basis functions
print "\n\n// Second derivatives of basis functions in interval (0,1)"
for N in range(0,max_order+1):
    for m in range(0,N+1):
        print("/** Second derivative of Lagrange basis polynomial of order %d which passes through the %d-th support point in interval (0,1). */" % (N,m))
        print("double basisFunctionSecondDerivative_%d_%d(const double s);"  % (N,m))
    print ""
