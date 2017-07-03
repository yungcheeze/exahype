#!/bin/bash
#
# Perform multicore speedup tests on Hamilton.
#
# Hamilton uses SLURM. SLURM supports array jobs.
#
# System specification(s):
#
# Hamilton 6 (x122 nodes):
#    2 x Intel Xeon E5-2650 v2 (Ivy Bridge) 8 cores, 2.6 GHz processors (16 cores per node)
#    64 GB DDR3 memory (4 GB per core)
#    the nodes are diskless
#    1 x TrueScale 4 x QDR single-port InfiniBand interconnect
#
# Hamilton 7 (x112 nodes):
#   2 x Intel Xeon E5-2650 v4 (Broadwell) 12 cores, 2.2 GHz processors (24 cores per node)
#   64 GB TruDDR4 memory
#   the nodes are diskless
#   1 x Intel OmniPath 100 Gb InfiniBand interconnect
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
 # T=(0.01 0.002 0.0005 0.0001) # p=3
 # if (( order == 5 )); then
 #   T=(0.003)                    # p=5
 # fi
 # if (( order == 7 )); then
 #   T=(0.001)                    # p=9
 # fi
 # if (( order == 9 )); then
 #   T=(0.0003)                   # p=9
 # fi


  # Create script
  script=multicore/hamilton.slurm-script
  newScript=multicore/hamilton-$prefix-p$order-n1-t1.slurm-script
  cp $script $newScript
 
  sed -i 's,'$project'-no-output-regular-0,'$prefix',g' $newScript

  sed -i 's,kernels=gen,kernels='$kernels',g' $newScript
  
  sed -i 's,p3,p'$order',g' $newScript
  sed -i 's,regular-0,'$mesh',g' $newScript

  sed -i 's,script=multicore/hamilton.slurm-script,script='$newScript',g' $newScript
  
  # Create spec files
  for coresPerTask in 24
  #for coresPerTask in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 48 # ham7
  #for coresPerTask in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 32 # ham6
  do
    spec=multicore/EulerADERDG-$io.exahype
    filename=multicore/$prefix-p$order-t1-c$coresPerTask
    newSpec=$filename'.exahype'

    cp $spec $newSpec

    if [[ "$kernels" == "opt" ]]; then
      sed -i -r "s,generic::fluxes::nonlinear,optimised::fluxes::nonlinear," $newSpec
    fi

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
