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
    # Gauss-Legendre nodes and weights.
    s, w = np.polynomial.legendre.leggauss(N+1)
    # Map onto (0,1).
    sGPN = 0.5*(s+1)
    wGPN = 0.5*w
    for m in range(0,N+1):
        s=sympy.symbols('s')
        print("double basisFunction_%d_%d(const double s) {" % (N,m))
        refphi = LagrangBasisPoly(s,N,m,sGPN.tolist())
        ret = "  return %s;" % sympy.simplify(refphi)
        ret = re.sub(r"s\*\*([0-9]+)", r"std::pow(s, \1)", ret)
        print ret
        print("}\n")

print "\n\n// First derivative of basis functions in interval (0,1)."
for N in range(0,max_order+1):
    # Gauss-Legendre nodes and weights.
    s, w = np.polynomial.legendre.leggauss(N+1)
    # Map onto (0,1).
    sGPN = 0.5*(s+1)
    wGPN = 0.5*w
    for m in range(0,N+1):
        s=sympy.symbols('s')
        refphi     = LagrangBasisPoly(s,N,m,sGPN.tolist())
        dds_refphi = sympy.diff(refphi, s)
        print("double basisFunctionFirstDerivative_%d_%d(const double s) {" % (N,m))
        ret = "  return %s;" % sympy.simplify(dds_refphi)
        ret = re.sub(r"s\*\*([0-9]+)", r"std::pow(s, \1)", ret)
        print ret
        print("}\n")

print "\n\n// Second derivative of basis functions in interval (0,1)."
for N in range(0,max_order+1):
    # Gauss-Legendre nodes and weights.
    s, w = np.polynomial.legendre.leggauss(N+1)
    # Map onto (0,1).
    sGPN = 0.5*(s+1)
    wGPN = 0.5*w
    for m in range(0,N+1):
        s=sympy.symbols('s')
        refphi       = LagrangBasisPoly(s,N,m,sGPN.tolist())
        d2ds2_refphi = sympy.simplify(sympy.diff(refphi, s, s))
        print("double basisFunctionSecondDerivative_%d_%d(const double s) {" % (N,m))
        ret = "  return %s;" % sympy.simplify(d2ds2_refphi)
        ret = re.sub(r"s\*\*([0-9]+)", r"std::pow(s, \1)", ret)
        print ret
        print("}\n")
