#!/usr/bin/env python
# See SRMHDBaseTest for more details

from SRMHDBaseTest import SRMHDTest, ConvergenceFrontend

test = SRMHDTest("ShockTubeSRMHDTest")

test.setExtend(width=1.0, height=1.0)
test.setInitialData("shocktube")

test.settings['ExaBinary'] = "../../../ApplicationExamples/SRMHD_LimitingADERDG/ExaHyPE-MHDSolver"

test.settings['ExaLimiter'] = 'Assumed to be enabled'


if __name__ == "__main__":
	ConvergenceFrontend(test, description=__doc__)

