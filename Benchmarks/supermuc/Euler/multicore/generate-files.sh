#!/bin/bash
#
# Perform multicore speedup tests on Hamilton.
#
# SuperMUC Phase 2 uses LOADLEVELER. 
# LOADLEVELER does not support array jobs.
#
# System specification:
#
# SuperMUC Phase 2 (x3072 nodes):
#   2x Haswell Xeon Processor E5-2697 v3 (2x14 cores)
#   64 GB memory
#   the nodes are diskless
#   Infinityband FDR14 interconnects
project=Euler

skipReductionInBatchedTimeSteps=on
batchFactor=0.8
hMax=( 0.03704 0.01235 0.00412 0.00138 0.00046 ) # 1/3^l ceiled with significance 1e-5
io=no-output # or output

kernels=gen # this is just an identifier; actual kernels must be chosen before building the executables

# Derived options

for fuseAlgorithmicSteps in "on" "off"
do
i=0
mesh=regular-$i
h=${hMax[i]}

prefix=$project-$io-$kernels
if [ "$fuseAlgorithmicSteps" == "on" ]; then
  prefix+="-fused"
else
  prefix+="-nonfused"
fi
prefix+="-$mesh"

for order in 3 5 7 9
do
  # SIMULATION END TIME
  T=( 0.01 0.00334 0.00112 0.00038 0.00013 )            # p=3
  if (( order == 5 )); then
    T=( 0.006364 0.002126 0.000713 0.000242 0.000083 )  # p=5; (2*3+1)/(2*order+1)*T_3 ceiled with sig. 1e-6
  fi
  if (( order == 7 )); then
    T=( 0.004667 0.001559 0.000523 0.000178 0.000061 )  # p=7
  fi
  if (( order == 9 )); then
    T=( 0.003685 0.001231 0.000413 0.00014 0.000048 )   # p=9
  fi
  t=${T[i]}

  for coresPerTask in 1 14 28
  do
    # Create script
    script=multicore/supermuc.load-leveler
    newScript=multicore/supermuc-$prefix-p$order-n1-t1-c$coresPerTask.load-leveler
    cp $script $newScript
 
    sed -i 's,coresPerTask=1,coresPerTask='$coresPerTask',g' $newScript
    
    sed -i 's,'$project'-no-output-regular-0,'$prefix',g' $newScript
    sed -i 's,p3,p'$order',g' $newScript
    sed -i 's,script=multicore/supermuc.load-leveler,script='$newScript',g' $newScript 

    # Create spec files  
    spec=multicore/$project-$io.exahype
    filename=multicore/$prefix-p$order-t1-c$coresPerTask 
    newSpec=$filename'.exahype'

    cp $spec $newSpec
    
    sed -i -r 's,end-time(\s*)=(\s*)(([0-9]|\.)*),end-time\1=\2'$t',' $newSpec
    sed -i -r 's,ranks_per_node:([0-9]+),ranks_per_node:1,' $newSpec 
    sed -i -r 's,cores(\s+)=(\s+)([0-9]+),cores\1=\2'$coresPerTask',' $newSpec
   
    sed -i -r 's,skip-reduction-in-batched-time-steps(\s*)=(\s*)(\w+),skip-reduction-in-batched-time-steps\1=\2'$skipReductionInBatchedTimeSteps',' $newSpec
    sed -i -r 's,timestep-batch-factor(\s*)=(\s*)(([0-9]|\.)+),timestep-batch-factor\1=\2'$batchFactor',' $newSpec
    sed -i -r 's,fuse-algorithmic-steps(\s*)=(\s*)(\w+),fuse-algorithmic-steps\1=\2'$fuseAlgorithmicSteps',' $newSpec
  
    sed -i -r 's,order(\s+)const(\s+)=(\s+)([0-9]+),order\1const\2=\3'$order',' $newSpec
    sed -i -r 's,maximum-mesh-size(\s*)=(\s*)(([0-9]|\.)*),maximum-mesh-size\1=\2'$h',' $newSpec
  done
done
done
