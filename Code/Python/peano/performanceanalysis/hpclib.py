"""
.. module:: hpclib
  :platform: Unix, Windows
  :synopsis: Contains functions to compute HPC measures.
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>

:synopsis: Contains functions to compute HPC measures.
"""

def compute_speedup(t,tref):
    """
    Args:
       t  (float[]): Array holding times.
       tref (float): Reference time.
    Returns:
       S (float[]) Speedup with respect to 'tref'.
    """
    n_t = len(t)
    S    = [0.0]*n_t
    
    for i in range(0,n_t):
        S[i] = tref/t[i]
    return S

def compute_speedup_2(t,tref):
    """
    Args:
       t (float[][]): Two-dimensional array holding times.
       tref (float): Reference time.
    Returns:
       S (float[][]) Speedup with respect to 'tref'.
    """
    n_ti = len(t)
    S  = [[]]*n_ti;
    for i in range(0,n_ti):
        S[i] = compute_speedup(t[i],tref)
    return S
