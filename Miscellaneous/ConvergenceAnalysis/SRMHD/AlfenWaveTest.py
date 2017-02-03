#!/usr/bin/env python
# See SRMHDBaseTest for more details

"""
SRMHD AlfenWave test.
"""

from SRMHDBaseTest import SRMHDTest, ConvergenceApplication

test = SRMHDTest("AlfenWaveTest")

test.setExtend(width=1.0, height=1.0)
test.setInitialData("AlfenWave")

# define a smaller test
until = test.until

#test.adaptiveRunRange = {
#	2: until(81), 3: until(0), 4: until(0),
#	5: until(0), 6: until(0), 7: until(0),
#	8: until(0), 9: until(0)
##}

# end  of smaller test

app = ConvergenceApplication(test, description=__doc__)


if __name__ == "__main__":
	app.parse_args()

