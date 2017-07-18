#!/usr/bin/python
# Python(2) script for computing post shock fluid
# pressure and velocity
# reference 1: https://en.wikipedia.org/wiki/Sod_shock_tube
# reference 2: https://gitlab.com/fantaz/Riemann_exact/tree/master
# reference 3: Toro, E., Riemann Solvers and Numerical Methods for Fluid Dynamics (4.3.3 Numerical tests)
from sympy import *
from mpmath import mp

mp.dps = 50

# Initial data
gamma=1.4

x_0 = 0.5;

rho_5 = 0.125; # right initial states
P_5   = 0.1;
u_5   = 0.0;

rho_1 = 1.0;    # left initial states
P_1   = 1.0;
u_1   = 0.0;

Gamma =(gamma-1)/(gamma+1)
Gamma2=Gamma**2

beta = (gamma-1)/(2*gamma)

print "constexpr double gamma     =%1.20f;"   % gamma
print "constexpr double Gamma     =%1.20f;"   % Gamma
print "constexpr double x_0       =%1.20f;\n" % x_0
print "constexpr double rho_5     =%1.20f;"   % rho_5
print "constexpr double P_5       =%1.20f;"   % P_5
print "constexpr double u_5       =%1.20f;"   % u_5
print "constexpr double rho_1     =%1.20f;"   % rho_1
print "constexpr double P_1       =%1.20f;"   % P_1
print "constexpr double u_1       =%1.20f;\n" % u_1

cs_1 = sqrt( gamma*P_1/rho_1 )
cs_5 = sqrt( gamma*P_5/rho_5 )

print "constexpr double cs_1       =%1.20f;"   % cs_1
print "constexpr double cs_5       =%1.20f;\n" % cs_5

P = symbols('P')
expressionForP_3 =  ( 
  (pow(P_1,beta)-pow(P,beta)) * sqrt( (1-Gamma2)*pow(P_1,gamma) / (Gamma2*rho_1) )
  -
  (P-P_5) * sqrt( (1-Gamma) / ( rho_5*(P + Gamma*P_5) ) )
)

lambdifiedExpressionForP_3 = lambdify(P, expressionForP_3, 'mpmath') 
P_3 = mpmath.findroot(lambdifiedExpressionForP_3, 0.3031)
P_4 = P_3

u_3   = u_5 + (P_3-P_5) / sqrt( rho_5/2 * ( (gamma+1)*P_3 + (gamma-1)*P_5 ) )
u_4   = u_3

u_2   =  (pow(P_1,beta)-pow(P_3,beta)) * sqrt( (1-Gamma2)*pow(P_1,gamma) / (Gamma2*rho_1) )

rho_3 = rho_1 * pow(P_3/P_1,1.0/gamma)
rho_4 = rho_5 * (P_4+Gamma*P_5)/(P_5+Gamma*P_4)

print "constexpr double u_2       =%1.20f;\n"   % u_2

cs_3 = sqrt( gamma*P_3/rho_3 )

print "constexpr double rho_3     =%1.20f;"   % rho_3
print "constexpr double P_3       =%1.20f;"   % P_3
print "constexpr double u_3       =%1.20f;"   % u_3
print "constexpr double cs_3      =%1.20f;\n" % cs_3

print "constexpr double rho_4     =%1.20f;"   % rho_4
print "constexpr double P_4       =%1.20f;"   % P_4
print "constexpr double u_4       =%1.20f;\n"   % u_4

u_shock = cs_5 * sqrt( 1+(gamma+1)/(2*gamma)*(P_4/P_5-1) )
print "constexpr double u_shock   =%1.20f;\n"   % u_shock
