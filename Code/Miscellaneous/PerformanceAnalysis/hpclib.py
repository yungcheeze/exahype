"""
.. module:: hpclib
  :platform: Unix, Windows, Mac
  :synopsis: Contains functions to compute HPC measures.
.. moduleauthor:: Dominic Etienne Charrier <dominic.e.charrier@durham.ac.uk>

:synopsis: Contains functions to compute HPC measures.
"""

def compute_speedup(t,t_ref):
    """
    Args:
       t  (float[]):
          Array holding times.
       t_ref (float):
          Reference time.
    Returns:
       Speedup (float[]) with respect to 't_ref'.
    """
    n_t = len(t)
    S    = [0.0]*n_t
    
    for i in range(0,n_t):
        S[i] = t_ref/t[i]
    return S

def compute_speedup_2(t,t_ref):
    """
    Args:
       t (float[][]):
          Two-dimensional array holding times.
       t_ref (float):
          Reference time.
    Returns:
       Speedup (float[][]) with respect to 't_ref'.
    """
    n_ti = len(t)
    S  = [[]]*n_ti;
    for i in range(0,n_ti):
        S[i] = compute_speedup(t[i],t_ref)
    return S
