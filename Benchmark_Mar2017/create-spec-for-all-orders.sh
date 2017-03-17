#!/bin/bash

hmaxs=(0.05 0.01 0.005 0.001)

for arrayid in 0 1 2 3
do

hmax=${hmaxs[arrayid]}

echo $hmax

for specfile in 'EulerFlow-output' 'EulerFlow-no-output'
do
for p in 3 4 5 6 7 8 9
do

newspecfile=$specfile'-regular-'$arrayid'-p'$p

cp $specfile'.exahype' $newspecfile'.exahype'
sed -i -r "s,order(\s*)const(\s*)=(\s*)([0-9]),order\1const\2=\3${p}," $newspecfile'.exahype'
sed -i -r "s,maximum-mesh-size(\s*)=(\s*)(([0-9]|\.)*),maximum-mesh-size\1=\2${hmax}," $newspecfile'.exahype'

done
done
done

