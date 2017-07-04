#!/usr/bin/python
# Python(2) script for computing post shock fluid
# pressure and velocity
#reference: http://www.phys.lsu.edu/~tohline/PHYS7412/sod.html
from sympy import *
from mpmath import mp

mp.dps = 50

P = symbols('P')

# Initial data
g=1.4
Pright=0.1
rright=0.125

m2=(g-1)/(g+1)

expressionForP =  ( 
  2*sqrt(g)/(g-1) * ( 1 - pow(P,(g-1)/(2*g)) ) 
  -
  ( P - Pright ) * sqrt( pow(1-m2,2) / ( rright * (P + m2*Pright ) ) )
)

lambdifiedExpressionForP = lambdify(P, expressionForP, 'mpmath') 
Proot = mpmath.findroot(lambdifiedExpressionForP, 0.3190)

vpost = 2*sqrt(g)/(g - 1) * ( 1 - pow(Proot,(g-1)/(2*g)) )

print "Ppost=%1.20f" % Proot
print "vpost=%1.20f" % vpost
