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
print "// Basis functions in interval (0,1)."
for N in range(0,max_order+1):
    for m in range(0,N+1):
        print("basisFunctions[%d][%d] = (UnivariateFunction) basisFunction_%d_%d;" % (N,m,N,m))

print "\n// First derivatives of basis functions in interval (0,1)."
for N in range(0,max_order+1):
    for m in range(0,N+1):
        print("basisFunctionFirstDerivatives[%d][%d] = (UnivariateFunction) basisFunctionFirstDerivative_%d_%d;" % (N,m,N,m))

print "\n// Second derivatives of basis functions in interval (0,1)."
for N in range(0,max_order+1):
    for m in range(0,N+1):
        print("basisFunctionSecondDerivatives[%d][%d] = (UnivariateFunction) basisFunctionSecondDerivative_%d_%d;" % (N,m,N,m))
