#!/bin/bash
hMax=(0.05 0.01 0.005 0.001)
times=(0.01 0.002 0.0005 0.0001)

i=0
mesh=regular-$i
h=${hMax[i]}
t=${times[i]}

order=5

sharedMem=None

skipReductionInBatchedTimeSteps=on
batchFactor=0.8

for io in 'output' 'no-output'
do
for coresPerTask in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 56
do 
  # Create script
  script=supermuc.load-leveler
  newScript=supermuc-$io-p$order-n1-t1-c$coresPerTask-$sharedMem.load-leveler
  cp $script $newScript
 
  sed -i -r 's,sharedMem=None,sharedMem='$sharedMem',' $newScript
  sed -i 's,EulerFlow_LimitingADERDG-no-output,EulerFlow_LimitingADERDG-'$io',g' $newScript
  sed -i 's,coresPerTask=1,coresPerTask='$coresPerTask',g' $newScript

  sed -i 's,p3,p'$order',g' $newScript
  sed -i 's,regular-0,'$mesh',g' $newScript

  # Create spec file
  spec=EulerFlow_LimitingADERDG-$io.exahype
  prefix=EulerFlow_LimitingADERDG-$io-p$order-$mesh-t1-c$coresPerTask
  newSpec=$prefix'.exahype'
  cp $spec $newSpec

  sed -i -r 's,end-time(\s*)=(\s*)(([0-9]|\.)*),end-time\1=\2'$t',' $newSpec
  sed -i -r 's,ranks_per_node:([0-9]+),ranks_per_node:1,g' $newSpec 
  sed -i -r 's,cores(\s+)=(\s+)([0-9]+),cores\1=\2'$coresPerTask',g' $newSpec
 
  sed -i -r 's,skip-reduction-in-batched-time-steps(\s*)=(\s*)(\w+),skip-reduction-in-batched-time-steps\1=\2'$skipReductionInBatchedTimeSteps',g' $newSpec
  sed -i -r 's,timestep-batch-factor(\s*)=(\s*)(([0-9]|\.)+),timestep-batch-factor\1=\2'$batchFactor',g' $newSpec
 
  sed -i -r 's,order(\s+)const(\s+)=(\s+)([0-9]+),order\1const\2=\3'$order',g' $newSpec
  sed -i -r 's,maximum-mesh-size(\s*)=(\s*)(([0-9]|\.)*),maximum-mesh-size\1=\2'$h',g' $newSpec

done
done
