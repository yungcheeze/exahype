#!/bin/bash

hmaxs=(0.05 0.01 0.005 0.001)
ts=(0.1 0.02 0.004 0.0008)

for i in 0 1 2 3
do

hmax=${hmaxs[i]}
t=${ts[i]}

echo $hmax


for specfile in 'EulerFlow-output' 'EulerFlow-no-output'
do
for p in 3 4 5 6 7 8 9
do

newspecfile=$specfile'-regular-'$i'-p'$p

cp $specfile'.exahype' $newspecfile'.exahype'

sed -i -r "s,order(\s*)const(\s*)=(\s*)([0-9]),order\1const\2=\3${p}," $newspecfile'.exahype'
sed -i -r "s,maximum-mesh-size(\s*)=(\s*)(([0-9]|\.)*),maximum-mesh-size\1=\2${hmax}," $newspecfile'.exahype'
sed -i -r "s,end-time(\s*)=(\s*)(([0-9]|\.)*),end-time\1=\2${t}," $newspecfile'.exahype'

done
done
done

