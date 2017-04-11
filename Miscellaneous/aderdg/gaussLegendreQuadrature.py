#!/usr/bin/python/

from numpy import *
 
def BaseFunc1D(xi, xin, N):
    """
    Computes the ADER-DG basis functions and their first derivative.
    
    Args:
       xi:
          The reference element point the basis functions are evaluated at.
          Here, xi refers to the greek letter that is often used as a reference element coordinate.
       xin:
          The reference element nodes corresponding to the nodal basis functions.
       N:
          Order of approximation corresponding to N+1 nodal basis functions.
    Returns:
       phi:
          Basis function values.
       phi_xi:
          First derivatives of the basis functions.
    """
    phi    = [1.]*(N+1) 
    phi_xi = [0.]*(N+1)
    for m in range(0,N+1):
        for j in range(0,N+1):
            if j == m:
                continue 
            phi[m] = phi[m]*(xi-xin[j])/(xin[m]-xin[j])
        for i in range(0,N+1):
            if i == m:
                continue
            tmp = 1.;
            for j in range(0,N+1):
                if j == i:
                    continue
                if j == m:
                    continue
                tmp = tmp*(xi-xin[j])/(xin[m]-xin[j])
            phi_xi[m] += tmp/(xin[m]-xin[i])
    return phi, phi_xi    

##################################################################
# Recursive generation of the Legendre polynomial of nodes n
def Legendre(n,x):
    x=array(x)
    if (n==0):
        return x*0+1.0
    elif (n==1):
        return x
    else:
        return ((2.0*n-1.0)*x*Legendre(n-1,x)-(n-1)*Legendre(n-2,x))/n
 
##################################################################
# Derivative of the Legendre polynomials
def DLegendre(n,x):
    x=array(x)
    if (n==0):
        return x*0
    elif (n==1):
        return x*0+1.0
    else:
        return (n/(x**2-1.0))*(x*Legendre(n,x)-Legendre(n-1,x))
##################################################################
# Roots of the polynomial obtained using Newton-Raphson method
def LegendreRoots(polynodes,tolerance=1e-20):
    if polynodes<2:
        err=1 # bad polynodes no roots can be found
    else:
        roots=[]
        # The polynomials are alternately even and odd functions. So we evaluate only half the number of roots. 
        for i in range(1,int(polynodes)/2 +1):
            x=cos(pi*(i-0.25)/(polynodes+0.5))
            error=10*tolerance
            iters=0
            while (error>tolerance) and (iters<1000):
                dx=-Legendre(polynodes,x)/DLegendre(polynodes,x)
                x=x+dx
                iters=iters+1
                error=abs(dx)
            roots.append(x)
        # Use symmetry to get the other roots
        roots=array(roots)
        if polynodes%2==0:
            roots=concatenate( (-1.0*roots, roots[::-1]) )
        else:
            roots=concatenate( (-1.0*roots, [0.0], roots[::-1]) )
        err=0 # successfully determined roots
    return [roots, err]
##################################################################
# Weight coefficients
def GaussLegendreWeights(polynodes):
    W=[]
    [xis,err]=LegendreRoots(polynodes)
    if err==0:
        W=2.0/( (1.0-xis**2)*(DLegendre(polynodes,xis)**2) )
        err=0
    else:
        err=1 # could not determine roots - so no weights
    return [W, xis, err]
##################################################################
# The integral value 
# func         : the integrand
# a, b         : lower and upper limits of the integral
# polynodes     : nodes of the Legendre polynomial to be used
#
def GaussLegendreQuadrature(func, polynodes, a, b):
    [Ws,xs, err]= GaussLegendreWeights(polynodes)
    if err==0:
        ans=(b-a)*0.5*sum( Ws*func( (b-a)*0.5*xs+ (b+a)*0.5 ) )
    else: 
        # (in case of error)
        err=1
        ans=None
    return [ans,err]
##################################################################
# The integrand - change as required
def func(x):
     nodes=2
     [xin,err]=LegendreRoots(nodes+1)
     
     phi_x,dphi_x = BaseFunc1D(x,xin,nodes)
     
     print phi_x
     
     return phi_x[0] * phi_x[0]
##################################################################
 
order=3
[Ws,xs,err]=GaussLegendreWeights(order+1)
if err==0:
    print "Order    : ", order
    print "Roots    : ", xs
    print "Weights  : ", Ws
else:
    print "Roots/Weights evaluation failed"

# Integrating the function
[ans,err]=GaussLegendreQuadrature(func , order+1, -1,1)
if err==0:
    print "Integral : ", ans
else:
    print "Integral evaluation failed"
