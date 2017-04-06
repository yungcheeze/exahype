#!/bin/bash
hMax=(0.5 0.1 0.05 0.01 0.005)

for order in 3 # 9
do

times=(0.1 0.05 0.01 0.005 0.001)
if (( order == 3 ));
then
  times=(1.0 0.5 0.1 0.05 0.01)
fi
# mesh
for i in 0 1 2 3 4
do
mesh=regular-$i
h=${hMax[i]}
t=${times[i]}

for fused in 'fused' 'nonfused'
do
for sharedMem in 'None' 'TBB'
do
for sharedMemIdentifier in 'dummy' 'autotuning' 'autotuning-without-learning'
do
for io in 'no-output' # output
do
#for coresPerTask in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 56
for coresPerTask in 1 2 4 7 14 28 56
do 
  # Create script
  script=supermuc.load-leveler
  newScript=supermuc-$io-$fused-$sharedMemIdentifier-p$order-n1-t1-c$coresPerTask-$sharedMem.load-leveler
  cp $script $newScript
 
  sed -i -r 's,sharedMem=None,sharedMem='$sharedMem',' $newScript
  sed -i 's,EulerFlow_LimitingADERDG-no-output,EulerFlow_LimitingADERDG-'$io'-'$fused'-'$sharedMemIdentifier',g' $newScript
  sed -i 's,coresPerTask=1,coresPerTask='$coresPerTask',g' $newScript 

  sed -i 's,p3,p'$order',g' $newScript
  sed -i 's,regular-0,'$mesh',g' $newScript
  
  # Create spec file
  spec=EulerFlow_LimitingADERDG-$io.exahype
  prefix=EulerFlow_LimitingADERDG-$io-$fused-$sharedMemIdentifier-p$order-$mesh-t1-c$coresPerTask
  newSpec=$prefix'.exahype'
  cp $spec $newSpec
  
  sed -i -r 's,end-time(\s*)=(\s*)(([0-9]|\.)*),end-time\1=\2'$t',' $newSpec 
 
  fuseSteps=off
  if [ '$fused'=='fused' ]; 
  then 
    fuseSteps=on   
  fi
  sed -i -r 's,fuse-algorithmic-steps(\s*)=(\s*)(\w+),fuse-algorithmic-steps\1=\2'$fuseSteps',g' $newSpec

  sed -i -r 's,identifier(\s+)=(\s+)(\w+),identifier\1=\2'$sharedMemIdentifier',g' $newSpec
  sed -i -r 's,cores(\s+)=(\s+)([0-9]+),cores\1=\2'$coresPerTask',g' $newSpec
  sed -i -r 's,properties-file(\s*)=(\s*)sharedmemory.'$h'.'$fused'.'$sharedMemIdentifier'.properties,g,' $newSpec
 
  sed -i -r 's,order(\s+)const(\s+)=(\s+)([0-9]+),order\1const\2=\3'$order',g' $newSpec
  sed -i -r 's,maximum-mesh-size(\s*)=(\s*)(([0-9]|\.)*),maximum-mesh-size\1=\2'$h',g' $newSpec

  echo 'created files '$newScript', '$newSpec

done
done
done
done 
done
done
done
