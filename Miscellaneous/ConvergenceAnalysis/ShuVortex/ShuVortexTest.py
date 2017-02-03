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
Run convergence tests: The EulerFlow ShuVortex.
Will setup environment variables to execute a templated specification file.
If called without arguments, it will run a series of runs in parallel.

Sample usages from the command line:

$ ./runShuVortex.py -p 2 -m 0.1
$ ./runShuVortex.py --all --wait
$ for p in 2 3 4; do ./runShuVortex.py -p $p -m 0.1 & done
"""

import sys, logging
logger = logging.getLogger("runShuVortex")
sys.path.append("../")

from libconvergence import PolyorderTest, ConvergenceApplication

test = PolyorderTest("EulerShuVortex")

# default values for ShuVortex simulations,
# cf. the paper of Michael Dumbser
test.polyorders = range(2,10)
test.width = 15.
test.depths = range(1,6)

test.computeCombinations()
until = test.until # shorthand

# big test. times are for single core runs, finest resolution:
test.adaptiveRunRange = {
	2: until(243), #  270min
	3: until(243), #  770min
	4: until(243), # 2000min (33h)
	5: until(81),  #  200min (243: > 2700min)
	6: until(81),  #  400min
	7: until(81),  #  800min
	8: until(27),  #   70min
	9: until(27)   #  130min
}


# smaller test, for analaysis of mistakes:
#test.adaptiveRunRange = { k: until(0) for k in range(2,10) }
#test.adaptiveRunRange[4] = until(81)

test.settings['ExaBinary'] = "../../../ApplicationExamples/EulerFlow/ExaHyPE-Euler"
# ONLY POLYORDER == 4
#test.settings['ExaBinaryTpl'] = test.settings['ExaBinary']

test.settings['ExaSpecfile'] = "ShuVortexConvergenceTpl.exahype"
# template to set up a queueing system
test.settings['QRUNTPL'] = "" # "mpirun -np 4" #"srun -n1 --partition=x-men --time=29:00:00 --mem=0 --job-name=p{ExapOrder}-m{ExaMeshSize}-ShuVortex"

# set initial data to use.
#settings['EXAHYPE_INITIALDATA']="MovingGauss2D"
test.settings['EXAHYPE_INITIALDATA']="ShuVortex"
# parameters for setting up the specfile
test.settings['ExaWidth'] = test.width
test.settings['ExaEndTime'] = 6.0
# parameters deciding how frequently output is made. As a first criterion,
# 1 output dump with the highest resolution is 250MB.
test.settings['ExaConvOutputRepeat'] = 0.1
test.settings['ExaVtkOutputRepeat'] = 0.2
# multi threads for the win
test.settings['ExaTbbCores'] = 1

# this is useful if you compiled with assertions
test.settings['EXAHYPE_SKIP_TESTS'] = True

app = ConvergenceApplication(test, description=__doc__)

if __name__ == "__main__":
	app.parse_args()

