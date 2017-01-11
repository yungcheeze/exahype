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
Run convergence tests: The SRMHD AlfenWave with Finite Volume.

Will setup environment variables to execute a templated specification file.
If called without arguments, it will run a series of runs in parallel.

Sample usages from the command line:

$ ./runThisTest.py -p 2 -m 0.1
$ ./runThisTest.py --all --wait
$ for p in 2 3 4; do ./runThisTest.py -p $p -m 0.1 & done
"""

import sys, logging
logger = logging.getLogger("runMhdFvAlfenTest")
sys.path.append("../libconvergence")

from convergence_test import PolyorderTest
from convergence_frontend import ConvergenceFrontend

test = PolyorderTest("MhdFvAlfenWaveTest")

# default values for ShuVortex simulations,
# cf. the paper of Michael Dumbser
test.polyorders = range(2,10)
test.width = 1.
test.depths = range(0,4)

test.computeCombinations()
until = test.until # shorthand

test.adaptiveRunRange = {
	2: until(27), 3: until(27), 4: until(27),
	5: until(27), 6: until(27), 7: until(27),
	8: until(27), 9: until(27)
}

test.settings['ExaBinary'] = "../../../ApplicationExamples/MHD_FV/ExaHyPE-MHDSolver"

test.settings['ExaSpecfile'] = "MHD_AlfenWaveConvergence.exahype"
# template to set up a queueing system
test.settings['QRUNTPL'] = "" #"srun -n1 --partition=x-men --time=29:00:00 --mem=0 --job-name=p{ExapOrder}-m{ExaMeshSize}-ShuVortex"

# Finite Volume:
# patch size of 14 cells in every direction
test.settings['PATCH_SIZE'] = 14

# set initial data to use.
#settings['EXAHYPE_INITIALDATA']="MovingGauss2D"
test.settings['EXAHYPE_INITIALDATA']="AlfenWaveThisIsNotUsed"
# parameters for setting up the specfile
test.settings['ExaWidth'] = test.width
test.settings['ExaEndTime'] = 10.0
# parameters deciding how frequently output is made. As a first criterion,
# 1 output dump with the highest resolution is 250MB.
test.settings['ExaConvOutputRepeat'] = 0.1
test.settings['ExaVtkOutputRepeat'] = 0.2
# single threaded in the moment.
test.settings['ExaTbbCores'] = 1

# this is useful if you compiled with assertions
test.settings['EXAHYPE_SKIP_TESTS'] = True


if __name__ == "__main__":
	ConvergenceFrontend(test, description=__doc__)

