#!/bin/bash
hMax=(0.05 0.01 0.005 0.001)
times=(0.01 0.002 0.0005 0.0001)

i=1
mesh=regular-$i
h=${hMax[i]}
t=${times[i]}

order=3
sharedMem=None

bufferSize=64
timeout=120

skipReductionInBatchedTimeSteps=on
batchFactor=0.8

for io in 'output' 'no-output'
do
for nodes in 10 28 82 
do
for tasksPerNode in 1 2 4 28
do  
  queue=micro
  if (( nodes > 20 )); then
    queue=general
  fi
  let tasks=$nodes*$tasksPerNode
  let coresPerTask=28/$tasksPerNode  

  # Create script
  script=supermuc.load-leveler
  newScript=supermuc-$io-p$order-n$nodes-t$tasksPerNode-c$coresPerTask.load-leveler
  cp $script $newScript
 
  sed -i -r 's,@ class(\s*)=(\s*)(\w+),@ class\1=\2'$queue',' $newScript
  sed -i -r 's,@ tasks_per_node(\s*)=(\s*)(([0-9]|\.)*),@ tasks_per_node\1=\2'$tasksPerNode',' $newScript
  sed -i -r 's,@ node(\s*)=(\s*)(([0-9]|\.)*),@ node\1=\2'$nodes',' $newScript
 
  sed -i 's,EulerFlow-no-output,EulerFlow-'$io',g' $newScript
  sed -i -r 's,sharedMem=None,sharedMem='$sharedMem',' $newScript

  sed -i 's,nodes=1,nodes='$nodes',g' $newScript
  sed -i 's,tasks=1,tasks='$tasks',g' $newScript
  sed -i 's,p3,p'$order',g' $newScript
  sed -i 's,regular-0,'$mesh',g' $newScript

  # Create spec file
  spec=EulerFlow-$io.exahype
  prefix=EulerFlow-$io-p$order-$mesh-t$tasksPerNode-c$coresPerTask
  newSpec=$prefix'.exahype'
  cp $spec $newSpec

  sed -i -r 's,end-time(\s*)=(\s*)(([0-9]|\.)*),end-time\1=\2'$t',' $newSpec
  sed -i -r 's,ranks_per_node:([0-9]+),ranks_per_node:'$tasksPerNode',g' $newSpec 
  sed -i -r 's,cores(\s+)=(\s+)([0-9]+),cores\1=\2'$coresPerTask',g' $newSpec
  sed -i -r 's,buffer-size(\s+)=(\s+)([0-9]+),buffer-size\1=\2'$bufferSize',g' $newSpec
  sed -i -r 's,timeout(\s+)=(\s+)([0-9]+),timeout\1=\2'$timeout',g' $newSpec
  
  sed -i -r 's,skip-reduction-in-batched-time-steps(\s*)=(\s*)(\w+),skip-reduction-in-batched-time-steps\1=\2'$skipReductionInBatchedTimeSteps',g' $newSpec
  sed -i -r 's,timestep-batch-factor(\s*)=(\s*)(([0-9]|\.)+),timestep-batch-factor\1=\2'$batchFactor',g' $newSpec

  sed -i -r 's,order(\s+)const(\s+)=(\s+)([0-9]+),order\1const\2=\3'$order',g' $newSpec
  sed -i -r 's,maximum-mesh-size(\s*)=(\s*)(([0-9]|\.)*),maximum-mesh-size\1=\2'$h',g' $newSpec

done
done
done
