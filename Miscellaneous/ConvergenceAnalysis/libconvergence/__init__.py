# The libconvergence package.
# (c) 2016,17 SvenK for ExaHyPE

"""
Libconvergence is a set of python classes to manage ExaHyPE runs
focussed on convergence studies. It also contains lot's of CLI frontend
code and analysis code for convergence studies.

It heavily uses pandas.
"""

# when from libconvergence import *:
#__all__ = ["echo", "surround", "reverse"]

# Shorthands for quicker use

from convergence_application import ConvergenceApplication
from convergence_table import ConvergenceReporter
from convergence_test import ConvergenceTest, PolyorderTest

# there is automatically present:
# libconvergence.convergence_application
# libconvergence.convergence_helpers      
# libconvergence.convergence_arghelpers
# libconvergence.convergence_table 
