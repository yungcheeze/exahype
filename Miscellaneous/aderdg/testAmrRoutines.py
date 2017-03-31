#!/usr/bin/python
"""
.. module:: main
  :platform: Unix, Windows, Mac
  :synopsis: Generates lookup table initialisation code up to a specific order.

:synopsis: Generates lookup table initialisation code up to a specific order.
"""
import numpy as np

import matplotlib.pylab
import matplotlib.pyplot as plt

from aderdg import *
from amrRoutines import *
from lagrangeInterpolation import *

from sympy.utilities.lambdify import lambdify

#-----------------------------------------------------------
# main
#-----------------------------------------------------------
maxOrder = 3;
minOrder = 3;
minDim   = 2;
maxDim   = 3;

order = minOrder
# Gauss-Legendre nodes and weights.
x, w = np.polynomial.legendre.leggauss(order+1)
# Map onto [0,1]
xGPN = 0.5*(x+1)
wGPN = 0.5*w
print(xGPN.tolist())

# Fine grid projectors.
fineGridProjector1d1 = assembleFineGridProjector1d(xGPN, 0, order)
fineGridProjector1d2 = assembleFineGridProjector1d(xGPN, 1, order)
fineGridProjector1d3 = assembleFineGridProjector1d(xGPN, 2, order)

factor = 9;
scaledXGPN = factor * xGPN

# sympy
x=sympy.symbols('x')

# subintervals
# 9*(0,1/3)
xGPN1 = 1./3. * scaledXGPN
Ly1 = [0.1, 0.25, 0.1, 0.25]
poly1 = LagrangePoly(x,xGPN1.tolist(),Ly1);
print("\nFine level Poly 1:")
print(Ly1)
print(sympy.simplify(poly1))
# 9 * (1/3,2/3)
xGPN2 = 1./3. * (scaledXGPN+1.*factor)
Ly2 = [0.15, 0.2, 0.1, 0.15]
poly2 = LagrangePoly(x,xGPN2.tolist(),Ly2);
print("\nFine level Poly 2:")
print(Ly2)
print(sympy.simplify(poly2))
# 9 * (2/3,1)
xGPN3 = 1./3. * (scaledXGPN+2.*factor) 
Ly3 = [0.15, 0.1,0.2, 0.05]
poly3 = LagrangePoly(x,xGPN3.tolist(),Ly3);
print("\nFine level Poly 3:")
print(Ly3)
print(sympy.simplify(poly3))

# perform restriction
coarseDof = np.zeros(order+1)
ones = np.ones(order+1)
singleLevelRestriction1d(coarseDof,np.array(Ly1),fineGridProjector1d1,wGPN,1,order)
singleLevelRestriction1d(coarseDof,np.array(Ly2),fineGridProjector1d2,wGPN,1,order)
singleLevelRestriction1d(coarseDof,np.array(Ly3),fineGridProjector1d3,wGPN,1,order)
print("\nCoarse DoF:")
print(coarseDof)

# interval
coarsePoly = LagrangePoly(x,scaledXGPN.tolist(),coarseDof.tolist());
print("\nCoarse Poly:")
print(sympy.simplify(coarsePoly))

# lambdify: express as callable python function
lam_poly1      = lambdify(x, poly1)
lam_poly2      = lambdify(x, poly2)
lam_poly3      = lambdify(x, poly3)
lam_coarsePoly = lambdify(x, coarsePoly)

# plot
t1      = np.arange(0.0,factor/3,0.1)
t2      = np.arange(factor/3.,factor*2./3.,0.1)
t3      = np.arange(factor*2./3.,factor,0.1)
tCoarse = np.arange(0,factor,0.1)

plt.plot(xGPN1,Ly1,antialiased=False,ls='None',marker='*')
plt.plot(t1,lam_poly1(t1))
plt.plot(xGPN2,Ly2,antialiased=False,ls='None',marker='*')
plt.plot(t2,lam_poly2(t2))
plt.plot(xGPN3,Ly3,antialiased=False,ls='None',marker='*')
plt.plot(t3,lam_poly3(t3))
#
plt.plot(scaledXGPN,coarseDof,antialiased=False,ls='None',marker='*')
plt.plot(tCoarse,lam_coarsePoly(tCoarse))
#plt.show()
