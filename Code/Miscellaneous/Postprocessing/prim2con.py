#!/usr/bin/python
#
# Convert a table of primitive variables to conserved ones. Got
# these kind of files from Olindo, looking like
#
#  x               p               rho            velx             anything
#  0.12500000E-02  0.10000000E+02  0.10000000E+01 -0.60000000E+00  0.15000000E+02
#  0.25000000E-02  0.10000000E+02  0.10000000E+01 -0.60000000E+00  0.15000000E+02
#
# I convert them to conserved quantities and write a similar table
# which is then suitable for comparison with exahype using vtkfinal1d.py.
# -- Sven

from pylab import *

prim = genfromtxt("solution-Mignone.dat", names=True, skip_header=2)

def prim2con(V):
	gamma = 5./3.
	Q = zeros_like(V)

	rho, vx, vy, vz, p = V
	v2 = vx**2 + vy**2 + vz**2
	if any(v2 > 1):
		print "Superluminal velocity in prim2con"
	lf = 1.0 / sqrt(1.0 - v2)
	gamma1 = gamma/(gamma-1.0)
	w      = rho + gamma1*p
	ww     = w*lf**2

	Q[0] = rho*lf
	Q[1] = ww*vx
	Q[2] = ww*vy
	Q[3] = ww*vz
	Q[4] = ww - p - Q[0]

	return Q

# conservative variables
p, velx, rho = [prim[c] for c in ('p','velx','rho')]
vely = zeros_like(velx)
velz = zeros_like(velx)

V = [rho,velx,vely,velz,p]
Q = prim2con(V)

cols = list(Q)
cols.insert(0, prim['x']) # spatial position
colnames = 'x d sx sy sz e'
savetxt("cons-Mignone.dat", transpose(cols), fmt='%.8e', header=colnames)
