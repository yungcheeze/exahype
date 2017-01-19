#!/usr/bin/env python
# See SRMHDBaseTest for more details

"""
SRMHD AlfenWave test.
"""

from SRMHDBaseTest import SRMHDTest, ConvergenceApplication

test = SRMHDTest("AlfenWaveTest")

test.setExtend(width=1.0, height=1.0)
test.setInitialData("AlfenWave")

app = ConvergenceApplication(test, description=__doc__)


if __name__ == "__main__":
	app.parse_args()

