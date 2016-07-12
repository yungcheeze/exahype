#!/usr/bin/python
"""
.. module:: main
  :platform: Unix, Windows, Mac
  :synopsis: Generates reference basis functions and their first and second derivatives up to a specific order.

:synopsis: Generates reference basis functions and their first and second derivatives up to a specific order.
"""
from lagrangeinterp import *
import numpy as np
import sympy
from sympy.printing import print_ccode
import re

order = 8
# Gauss-Legendre nodes and weights.
s, w = np.polynomial.legendre.leggauss(order+1)
# Map onto (0,1).
sGPN = 0.5*(s+1)
wGPN = 0.5*w

# Print the 1d reference basis functions#
print "################################################################################"
print "Basis functions in reference interval (0,1) for order %d:" % order
print "################################################################################"
for m in range(0,order+1):
    s=sympy.symbols('s')
    print("double basis_function_%d_%d (s) {" % (order,m))
    refphi = LagrangBasisPoly(s,order,m,sGPN.tolist())
    ret = "  return %s;" % sympy.simplify(refphi)
    ret = re.sub(r"s\*\*([0-9]+)", r"std::pow(s, \1)", ret)
    print ret
    print("}")

# Print the first derivative of the 1d reference basis functions
print "\n################################################################################"
print "1st derivative of basis functions for order %d:" % order
print "################################################################################"
for m in range(0,order+1):
    s=sympy.symbols('s')
    refphi     = LagrangBasisPoly(s,order,m,sGPN.tolist())
    dds_refphi = sympy.diff(refphi, s)
    print("double basis_first_derivative_%d_%d (s) {" % (order,m))
    ret = "  return %s;" % sympy.simplify(dds_refphi)
    ret = re.sub(r"s\*\*([0-9]+)", r"std::pow(s, \1)", ret)
    print ret
    print("}")

## Print the second derivative of the 1d reference basis functions
print "\n################################################################################"
print "2nd derivative of basis functions for order %d:" % order
print "################################################################################"
for m in range(0,order+1):
    s=sympy.symbols('s')
    refphi       = LagrangBasisPoly(s,order,m,sGPN.tolist())
    d2ds2_refphi = sympy.simplify(sympy.diff(refphi, s, s))
    print("double basis_second_derivative_%d_%d (s) {" % (order,m))
    ret = "  return %s;" % sympy.simplify(d2ds2_refphi)
    ret = re.sub(r"s\*\*([0-9]+)", r"std::pow(s, \1)", ret)
    print ret
    print("}")
