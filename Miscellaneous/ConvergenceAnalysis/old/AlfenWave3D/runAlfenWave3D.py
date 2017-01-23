#!/usr/bin/env python
# Run convergence tests.

"""
Run convergence tests: The MHD AlfenWave
See ShuVortex to read more.
"""

import sys, logging
logger = logging.getLogger("runShuVortex")
sys.path.append("../libconvergence")

from convergence_starter import PolyorderTest, convergenceFrontend

test = PolyorderTest("SRMHD3D-AlfenWave")

# 3D MHD Wave
# with default MHD application but in 3D
test.polyorders = range(2,10)
test.width = 1.0
test.depths = range(1,6)

# Actually there is a problem: The depths are way to large for a 3D simulation
# on a single core. This is why most of the simulations fail.
test.computeCombinations()

until = test.until # shorthand
test.adaptiveRunRange = {
	2: until(243), 3: until(243), 4: until(243),
	5: until(243), 6: until(81), 7: until(81),
	8: until(27), 9: until(27)
}

test.settings['ExaBinary'] = "../../../ApplicationExamples/MHD/ExaHyPE-MHD"
test.settings['ExaSpecfile'] = "MHD_AlfenWave3DConvergence.exahype"

# template to set up a queueing system
# Comment out if you don't want to use it.
#test.settings['QRUNTPL'] = "srun -n1 --partition=x-men --time=29:00:00 --mem=0 --job-name=MHD3D-p{ExapOrder}-m{ExaMeshSize}-ConvergenceStudies"


# set initial data to use.
test.settings['EXAHYPE_INITIALDATA']="AlfenWave"
# parameters for setting up the specfile
test.settings['ExaWidth']=test.width
test.settings['ExaEndTime']=12.0 # MHD LONG RUN
# parameters deciding how frequently output is made. As a first criterion,
# 1 output dump with the highest resolution is 250MB.
test.settings['ExaConvOutputRepeat']=0.1
test.settings['ExaVtkOutputRepeat']=0.5
# single threaded in the moment.
test.settings['ExaTbbCores'] = 1

test.settings['ExaVtkFormat'] = "Legendre::vertices"

# this is useful if you compiled with assertions
test.settings['EXAHYPE_SKIP_TESTS'] = "True"

if __name__ == "__main__":
	convergenceFrontend(test, description=__doc__)


