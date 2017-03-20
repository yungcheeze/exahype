#!/bin/bash

hmaxs=(0.05 0.01 0.005 0.001)
ts=(0.01 0.002 0.0001 0.00002)

for i in 0 1 2 3
do

hmax=${hmaxs[i]}
t=${ts[i]}

echo $hmax


for specfile in $1 $2
do
for p in 3 7 9
do

newspecfile=$specfile'-regular-'$i'-p'$p

cp $specfile'.exahype' $newspecfile'.exahype'

sed -i -r "s,order(\s*)const(\s*)=(\s*)([0-9]),order\1const\2=\3${p}," $newspecfile'.exahype'
sed -i -r "s,maximum-mesh-size(\s*)=(\s*)(([0-9]|\.)*),maximum-mesh-size\1=\2${hmax}," $newspecfile'.exahype'
sed -i -r "s,end-time(\s*)=(\s*)(([0-9]|\.)*),end-time\1=\2${t}," $newspecfile'.exahype'

done
done
done

