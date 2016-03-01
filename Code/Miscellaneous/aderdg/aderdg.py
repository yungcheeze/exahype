import numpy as np
from numpy import linalg

"""
.. module:: aderdg
  :platform: Unix, Windows, Mac
  :synopsis: Provides routines to compute ADER-DG basis functions and operators on the unit cube.
.. moduleauthor:: Angelika Schwarz <angelika.schwarz@tum.de>

:synopsis: Provides routines to compute ADER-DG basis functions and operators on the unit cube.
"""

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


def assembleStiffnessMatrix(xGPN, wGPN, N):
    """
    Computes the (reference) element stiffness matrix for an approximation of
    order N.

    Args:
       xGPN:
          Gauss-Legendre nodes (N nodes).
       wGPN:
          N Gauss-Legendre weights  (N weights).
       N:
          Order of approximation corresponding to N+1 nodal basis functions.
    Returns:
       K_xi:
          The (reference) element stiffness matrix.
    """
    # init matrix with zero
    Kxi = [[0 for _ in range(N+1)] for _ in range(N+1)]
     
    for i in range(0,N+1):
        phi, phi_xi = BaseFunc1D(xGPN[i], xGPN, N)
        for k in range(0,N+1):
            for l in range(0,N+1):
                Kxi[k][l] += wGPN[i]*phi_xi[k]*phi[l] 
        
    return Kxi

def assembleMassMatrix(xGPN, wGPN, N):
    """
    Computes the (reference) element mass matrix for an approximation of
    order N.

    Args:
       xGPN:
          Gauss-Legendre nodes (N nodes).
       wGPN:
          N Gauss-Legendre weights (N weights).
       N:
          Order of approximation corresponding to N+1 nodal basis functions.
    Returns:
       M_xi:
          The (reference) element mass matrix.
    """
    # init matrix with zeros
    MM = [[0 for _ in range(N+1)] for _ in range(N+1)]
    
    for i in range(0,N+1):
        phi, _ = BaseFunc1D(xGPN[i], xGPN, N)
        for k in range(0,N+1):
            for l in range(0,N+1):
                MM[k][l] += wGPN[i]*phi[k]*phi[l]
      
    return MM

def assembleK1(Kxi, xGPN, N):
    """
    Computes the difference between the reference element mass operator 
    evaluated at point xi=1.0 and the element stiffness matrix.
    
    Args:
       K_xi:
          The (reference) element stiffness matrix for a approximation of 
          order N.
       xGPN:
          N Gauss-Legendre nodes (N nodes).
       N:
          Order of approximation corresponding to N+1 nodal basis functions.
    Returns:
       K1:
          <unknown>
    """
    phi1, _ = BaseFunc1D(1.0, xGPN, N)
    FRm = [[0 for _ in range(N+1)] for _ in range(N+1)]
    
    for k in range(0, N+1):
        for l in range(0, N+1):
            FRm[k][l] = phi1[k]*phi1[l] 
    
    K1 = np.subtract(FRm,Kxi)
    return K1        
        
def assembleTimeFluxMatrixF0(xGPN, N):
    """
    Evaluates the reference basis functions at point xi=0.0.

    Args:
       xGPN:
          Gauss-Legendre nodes (N nodes).
       N:
          Order of approximation corresponding to N+1 nodal basis functions.
    Returns:
       phi:
          The reference basis functions evaluated at point xi=0.0.
    """
    phi, _ = BaseFunc1D(0.0, xGPN, N)
    return phi         
        
def assembleDiscreteDerivativeOperator(MM,Kxi):
    """
    Computes some derivative values for debugging purposes.

    Args:
       MM:
          The (reference) element mass matrix for a approximation of 
          order N.
       Kxi:
          The (reference) element stiffness matrix for a approximation of 
          order N.
       
    Returns:
       dudx:
          Derivative values for debugging purposes.
    """
    dudx = np.dot(linalg.inv(MM),np.transpose(Kxi))
    return dudx
    
def assembleSubOutputMatrix(xGPN, N, dim):
    """
    Transforms the degrees of freedom located at the non-equidistant Gauss-Legendre 
    to degrees of freedoms hold at nodes of an uniform grid over (0,1)^dim.
    
    Let us denote by P the 1D projection operator (1DsubOutputMatrix). The uniform DoF 
    are computed according to:
    
    2D:
       u^eq_ij =  sum_{ij;mn} P_im P_jn u^DG_mn
    3D:
       u^eq_ijk =  sum_{ijk;mnp} P_im P_jn P_kp u^DG_mnp
    
    Args:
       xGPN:
          Gauss-Legendre nodes (N nodes).
       N:
          Order of approximation corresponding to N+1 nodal basis functions.
       dim:
          Space dimension.
    Returns:
       subOutputMatrix:
          The correspondng degrees of freedom located at an uniform grid over (0,1)^dim.
    """
    subOutputMatrix = [[0 for _ in range((N+1)**dim)] for _ in range((N+1)**dim)] # 16 x 16
    subxi = np.linspace(0.0, 1.0, num=(N+1))
    cnt = 0 
    
    if dim == 3:
        for k in range(0, N+1):
            for j in range(0, N+1):
                for i in range(0, N+1):
                    phi_i, _ = BaseFunc1D(subxi[i], xGPN, N)
                    phi_j, _ = BaseFunc1D(subxi[j], xGPN, N)
                    phi_k, _ = BaseFunc1D(subxi[k], xGPN, N)
                    count = 0
                    for kk in range(0, N+1):
                        for jj in range(0, N+1):
                            for ii in range(0, N+1):
                                aux = [phi_i[ii],phi_j[jj],phi_k[kk]]
                                val = reduce(lambda x, y: x*y, aux)
                                subOutputMatrix[count][cnt] = val
    
                                
                                count +=1
                    cnt+=1
    elif dim == 2:
        for j in range(0, N+1): # eq basis
            for i in range(0, N+1):
                phi_i, _ = BaseFunc1D(subxi[i], xGPN, N)
                phi_j, _ = BaseFunc1D(subxi[j], xGPN, N)
                count = 0
                for jj in range(0, N+1): # DG basis
                    for ii in range(0, N+1):
                        aux = [phi_i[ii],phi_j[jj]]
                        val = reduce(lambda x, y: x*y, aux)
                        subOutputMatrix[count][cnt] = val # count belongs to dg basis
                        
                        count +=1
                cnt+=1  
    else:
        sys.stderr.write("Dimension not supported. Suboutputmatrix not generated.")          
    return subOutputMatrix

def assembleEquidistantGridProjector1D(xGPN, N, dim):
    """
    Transforms the degrees of freedom located at the non-equidistant Gauss-Legendre
    nodes to degrees of freedoms located at nodes of an equidistant grid over (0,1).
    
    Let us denote by P the 1D projection operator (=1D subOutputMatrix). The equidistant DoF 
    are computed according to:
    
    u^eq_i =  sum_{m} P_im u^DG_m
    
    Args:
       xGPN:
          Gauss-Legendre nodes (N nodes).
       N:
          Order of approximation corresponding to N+1 nodal basis functions.
       dim:
          Space dimension.
    Returns:
       subOutputMatrix:
          The corresponding degrees of freedom located at nodes of an equidistant grid over (0,1).
    """
    subOutputMatrix = [[0 for _ in range((N+1))] for _ in range((N+1))] # 4 x 4 for N=3
    subxi = np.linspace(0.0, 1.0, num=(N+1))
    
    for i in range(0, N+1): # Eq. basis
        phi_i, _ = BaseFunc1D(subxi[i], xGPN, N) # comma is important here
        for m in range(0, N+1): # DG basis
            subOutputMatrix[m][i] = phi_i[m]
    return subOutputMatrix

def assembleEquidistantGridProjector2D(xGPN, N, dim):
    """
    Transforms the degrees of freedom located at the non-equidistant Gauss-Legendre
    nodes to degrees of freedoms located at nodes of an equidistant grid over (0,1).
    
    Let us denote by P the 1D projection operator (=1D subOutputMatrix). The equidistant DoF 
    are computed according to:
    
    u^eq_ij =  sum_{ij;mn} P_im P_jn u^DG_mn
    
    Args:
       xGPN:
          Gauss-Legendre nodes (N nodes).
       N:
          Order of approximation corresponding to N+1 nodal basis functions.
       dim:
          Space dimension.
    Returns:
       subOutputMatrix:
          The corresponding degrees of freedom located at nodes of an equidistant grid over (0,1).
    """
    subOutputMatrix = [[0 for _ in range((N+1)**2)] for _ in range((N+1)**2)] # 16 x 16 for N=3
    
    projector1D = assembleEquidistantGridProjector1D(xGPN, N, dim)
    
    for n in range(0, N+1): # dg basis
        for m in range(0, N+1):
            dofIndexDG = m +n*(N+1)
            for j in range(0, N+1): # eq basis
                for i in range(0, N+1):
                    dofIndexEq = i + j*(N+1)
                    subOutputMatrix[dofIndexDG][dofIndexEq] = projector1D[n][j] * projector1D[m][i]
    return subOutputMatrix

def assembleEquidistantGridProjector3D(xGPN, N, dim):
    """
    Transforms the degrees of freedom located at the non-equidistant Gauss-Legendre
    nodes to degrees of freedoms located at nodes of an equidistant grid over (0,1).
    
    Let us denote by P the 1D projection operator (=1D subOutputMatrix). The equidistant DoF 
    are computed according to:
    
    u^eq_ijk =  sum_{ijk;mnp} P_im P_jn P_kp u^DG_mnp
    
    Args:
       xGPN:
          Gauss-Legendre nodes (N nodes).
       N:
          Order of approximation corresponding to N+1 nodal basis functions.
       dim:
          Space dimension.
    Returns:
       subOutputMatrix:
          The corresponding degrees of freedom located at nodes of an equidistant grid over (0,1).
    """
    subOutputMatrix = [[0 for _ in range((N+1)**3)] for _ in range((N+1)**3)]
    
    projector1D = assembleEquidistantGridProjector1D(xGPN, N, dim)
    
    for p in range(0, N+1):
        for n in range(0, N+1):
            for m in range(0, N+1):
                dofIndexDG = m +n*(N+1)+p*(N+1)**2
                for k in range(0, N+1):
                    for j in range(0, N+1):
                        for i in range(0, N+1):
                            dofIndexEq = i +j*(N+1)+k*(N+1)**2
                            subOutputMatrix[dofIndexDG][dofIndexEq] = projector1D[p][k] * projector1D[n][j] * projector1D[m][i]
    return subOutputMatrix
