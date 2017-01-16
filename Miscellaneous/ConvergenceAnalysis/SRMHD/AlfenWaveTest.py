#!/usr/bin/env python
# See SRMHDBaseTest for more details

from SRMHDBaseTest import SRMHDTest, ConvergenceFrontend

test = SRMHDTest("AlfenWaveTest")

test.setExtend(width=1.0, height=1.0)
test.setInitialData("AlfenWave")


if __name__ == "__main__":
	ConvergenceFrontend(test, description=__doc__)

