#!/usr/bin/env python
#
# Run convergence tests.
#
# This is basically python version agnostic and should run with both python2
# and python3.
#
# ExaHyPE 2016, SvenK
#

"""
Run convergence tests: The Optimized EulerFlow ShuVortex
See the regular ShuVortex for comparison.
"""

import sys, logging
logger = logging.getLogger("runOptimizedShuVortex")
sys.path.append("../libconvergence")

from convergence_starter import PolyorderTest, convergenceFrontend

test = PolyorderTest("OptimizedEulerShuVortex")

# default values for ShuVortex simulations,
# cf. the paper of Michael Dumbser
test.polyorders = range(2,10)
test.width = 2.
test.depths = range(1,6)

test.computeCombinations()
until = test.until # shorthand

test.adaptiveRunRange = {
	2: until(243), 3: until(243), 4: until(243),
	5: until(243), 6: until(81), 7: until(81),
	8: until(27), 9: until(27)
}

test.settings['ExaBinary'] = "../../../ApplicationExamples/OptimisedKernel_Euler/ExaHyPE-Euler"

test.settings['ExaSpecfile'] = "OptimizedShuVortexConvTpl.exahype"
# template to set up a queueing system
test.settings['QRUNTPL'] = "srun -n1 --partition=x-men --time=29:00:00 --mem=0 --job-name=p{ExapOrder}-m{ExaMeshSize}-OptimizedShuVortex"

# set initial data to use.
#settings['EXAHYPE_INITIALDATA']="MovingGauss2D"
test.settings['EXAHYPE_INITIALDATA']="ShuVortex"
# parameters for setting up the specfile
test.settings['ExaWidth'] = test.width
test.settings['ExaEndTime'] = 2.0
# parameters deciding how frequently output is made. As a first criterion,
# 1 output dump with the highest resolution is 250MB.
test.settings['ExaConvOutputRepeat'] = 0.1
test.settings['ExaVtkOutputRepeat'] = 0.5
# single threaded in the moment.
test.settings['ExaTbbCores'] = 1

# this is useful if you compiled with assertions
test.settings['EXAHYPE_SKIP_TESTS'] = True


if __name__ == "__main__":
	convergenceFrontend(test, description=__doc__)

