#!/usr/bin/python/

from numpy import *
 
def singleLevelRestriction1d(coarseDof,fineDof,fineGridProjector1d,weights,Nvar,N):
    """
    Restricts the fine level degrees of freedom (DoF) up to the coarse level
    and adds them to the coarse level DoF.
    
    Args:
       coarseDof: array
          Degrees of freedom on the coarse level.
       fineDof: array
          Degrees of freedom on the fine level.
       fineGridProjector1d: list of list of strings
          Operator that is used to project the fine
          grid degrees of freedom.
       weights:
          Quadrature weights.
       Nvar:
          number of variables
       N:
          Order of approximation.
    Returns:
       coarseDof:
          Degrees of freedom on the coarse level including the
          contributions of the fineDof.
    """
    for m1 in range(0,N+1):
        for ivar in range(0,Nvar):
            mNodeIndex     = m1;
            mDofStartIndex = mNodeIndex * Nvar;

            for n1 in range(0,N+1):
               nNodeIndex     = n1;
               nDofStartIndex = nNodeIndex * Nvar;
               
               coarseDof[mDofStartIndex+ivar] += float(weights[n1]) * float(fineGridProjector1d[m1][n1]) * fineDof[nDofStartIndex + ivar] / float(weights[m1]) / 3.;
    return coarseDof

def singleLevelProlongation1d(fineDof,coarseDof,fineGridProjector1d,weights,Nvar,N):
    """
    Prolongates the coarse level DoF down to the fine level.
    
    Args:
       fineDof:
          Degrees of freedom on the fine level.
       coarseDof:
          Degrees of freedom on the coarse level.
       fineGridProjector1d:
          Operators that are used to project the fine
          grid degrees of freedom.
       weights:
          Quadrature weights.
       N:
          Order of approximation.
    Returns:
       coarseDof:
          Degrees of freedom on the coarse level including the
          contributions of the fineDof.
    """
    for m1 in range(0,N+1):
        for ivar in range(0,Nvar):
            mNodeIndex     = m1;
            mDofStartIndex = mNodeIndex * Nvar;

            for n1 in range(0,N+1):
               nNodeIndex     = n1;
               nDofStartIndex = nNodeIndex * Nvar;
               
               fineDof[mDofStartIndex+ivar] += coarseDof[nDofStartIndex + ivar] * float(fineGridProjector1d[n1][m1]);
    return fineDof
