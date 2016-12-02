#!/bin/bash
#
# a script to extract in lines like
#
# 2016-09-27 17:42:51  88832.4      [beast],rank:0 info         exahype::mappings::LoadBalancing::endIteration(State)   memoryUsage    =4136 MB
# 2016-09-27 17:42:51  88892.6      [beast],rank:0 info         exahype::runners::Runner::startNewTimeStep(...)         step 1282 t_min          =0.297167
# 2016-09-27 17:42:51  88892.6      [beast],rank:0 info         exahype::runners::Runner::startNewTimeStep(...)           dt_min         =0.000231799
# 2016-09-27 17:42:51  88898.3      [beast],rank:0 info         exahype::mappings::LoadBalancing::endIteration(State)   memoryUsage    =4136 MB
#
# the dt_min value.
#

for sim in simulations/p*; do
	echo "$sim"
	cd $sim
	SIMLOG=$(ls *log)
	cat $SIMLOG | grep dt_min | awk '{ print $8; }' | sed 's/=//' > dtminlog.asc
	cd - &>/dev/null
done

# we can now feed the asc files into gnuplot.
