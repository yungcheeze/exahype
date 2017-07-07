#!/bin/bash
project=EulerADERDG

skipReductionInBatchedTimeSteps=on
batchFactor=0.8
hMax=(0.05 0.01 0.005 0.001)
T=(0.2 0.2 0.0005 0.0001) # p=3
io=no-output # or output

kernels=gen

# Derived options

for fuseAlgorithmicSteps in "on" "off"
do
  prefix=$project-$io

  i=0
  mesh=regular-$i
  h=${hMax[i]}
  t=${T[i]}

  prefix+=-$kernels-$mesh


  if [ "$fuseAlgorithmicSteps" == "on" ]; then
   prefix+="-fused"
  else
   prefix+="-nonfused"
  fi

  # prefix= project-io-kernels-mesh-i-
  # regex = .+-.+-(\w+)-.+-.+-

  for order in 3 5 7 9
  do
    T=(0.01 0.002 0.0005 0.0001)   # p=3
    if (( order == 5 )); then
      T=(0.003)                    # p=5
    fi
    if (( order == 7 )); then
      T=(0.001)                    # p=7
    fi
    if (( order == 9 )); then
      T=(0.0003)                   # p=9
    fi
    
    for coresPerTask in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 56
    do 
      # Create script
      script=multicore/supermuc.load-leveler
      newScript=multicore/supermuc-$prefix-p$order-n1-t1-c$coresPerTask.load-leveler
      cp $script $newScript
     
      sed -i -r 's,sharedMem=None,sharedMem='$sharedMem',' $newScript
      sed -i 's,EulerADERDG-no-output,EulerADERDG-'$io',g' $newScript
      sed -i 's,coresPerTask=1,coresPerTask='$coresPerTask',g' $newScript 

      sed -i 's,p3,p'$order',g' $newScript
      sed -i 's,regular-0,'$mesh',g' $newScript

      sed -i 's,script=supermuc.load-leveler,script='$newScript',g' $newScript

      # Create spec file
      spec=EulerADERDG-$io.exahype
      prefix=EulerADERDG-$io-p$order-$mesh-t1-c$coresPerTask
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
done
