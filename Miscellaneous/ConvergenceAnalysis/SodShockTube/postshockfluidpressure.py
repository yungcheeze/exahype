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

x0 = 0.5;

rho_r = 0.125;
P_r   = 0.1;
u_r   = 0;

rho_l = 1;
P_l   = 1;
u_l   = 0;

mu2=(g-1)/(g+1)
mu = sqrt( mu2 );

print "constexpr double gamma     =%1.20f;" % g
print "constexpr double mu        =%1.20f;" % mu
print "constexpr double x0        =%1.20f;\n" % x0
print "constexpr double rho_r     =%1.20f;" % rho_r
print "constexpr double P_r       =%1.20f;" % P_r
print "constexpr double u_r       =%1.20f;" % u_r
print "constexpr double rho_l     =%1.20f;" % rho_l
print "constexpr double P_l       =%1.20f;" % P_l
print "constexpr double u_l       =%1.20f;\n" % u_l

c_l = sqrt( g*P_l/rho_l )
c_r = sqrt( g*P_r/rho_r )

print "constexpr double c_l       =%1.20f;" % c_l
print "constexpr double c_r       =%1.20f;\n" % c_r

expressionForP =  ( 
  2*sqrt(g)/(g-1) * ( 1 - pow(P,(g-1)/(2*g)) ) 
  -
  ( P - P_r ) * sqrt( pow(1-mu2,2) / ( rho_r * (P + mu2*P_r ) ) )
)

lambdifiedExpressionForP = lambdify(P, expressionForP, 'mpmath') 
P_post = mpmath.findroot(lambdifiedExpressionForP, 0.3190)

v_post = 2*sqrt(g)/(g - 1) * ( 1 - pow(P_post,(g-1)/(2*g)) )

rho_l = 1;
P_l = 1;
u_l = 0;

rho_r = 0.125;
P_r = 0.1;
u_r = 0;


print "constexpr double P_post    =%1.20f;" % P_post
print "constexpr double v_post    =%1.20f;" % v_post

rho_l = 1;
P_l = 1;
u_l = 0;

rho_r = 0.125;
P_r = 0.1;
u_r = 0;

rho_post   = rho_r*(( (P_post/P_r) + mu2 )/(1 + mu2*(P_post/P_r)));
v_shock    = v_post*((rho_post/rho_r)/( (rho_post/rho_r) - 1.0));
rho_middle = (rho_l)*pow((P_post/P_l),1.0/g);

print "constexpr double rho_post  =%1.20f;" % rho_post
print "constexpr double v_shock   =%1.20f;" % v_shock
print "constexpr double rho_middle=%1.20f;" % rho_middle
