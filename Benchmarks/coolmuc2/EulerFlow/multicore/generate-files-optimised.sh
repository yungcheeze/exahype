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

for order in 3 5 7 9
do
  hMax=(0.05 0.01 0.005 0.001)
  if [ "$order" -le 3 ];  then
    times=(0.01 0.002 0.0005 0.0001)
  elif [ "$order" -le 5 ];  then
    times=(0.003)
  elif [ "$order" -le 7 ];  then
    times=(0.001)
  else
    times=(0.0006) 
  fi

  i=0
  mesh=regular-$i
  h=${hMax[i]}
  t=${times[i]}


  skipReductionInBatchedTimeSteps=on
  batchFactor=0.8

  #for io in 'output' 'no-output'
  for io in 'no-output'
  do
    #for coresPerTask in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 48 # ham7
    #for coresPerTask in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 32 # ham6
    for coresPerTask in 1 2 4 14 28 56
    do
      # Create script
      script=coolmuc2.slurm-script
      newScript=coolmuc2-$io-p$order-n1-t1-c$coresPerTask-opt.slurm-script
      cp $script $newScript
     
      sed -i 's,opt=gen,opt=opt,g' $newScript
      sed -i 's,coresPerTask=1,coresPerTask='$coresPerTask',g' $newScript
      sed -i 's,EulerFlow-no-output,EulerFlow-'$io',g' $newScript

      sed -i 's,p3,p'$order',g' $newScript
      sed -i 's,regular-0,'$mesh',g' $newScript

      sed -i 's,script=coolmuc2.slurm-script,script='$newScript',g' $newScript  
      # Create spec files
      spec=EulerFlow-$io.exahype
      prefix=EulerFlow-$io-p$order-$mesh-t1-c$coresPerTask-opt
      newSpec=$prefix'.exahype'
      cp $spec $newSpec
      
      sed -i -r 's,generic::fluxes::nonlinear,optimised::fluxes::nonlinear,' $newSpec
      sed -i -r 's,architecture(\s+)const(\s+)=(\s+)(\w+),architecture\1const\2=\3hsw,g' $newSpec

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
