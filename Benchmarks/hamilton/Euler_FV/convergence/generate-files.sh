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
project=Euler_FV

skipReductionInBatchedTimeSteps=on
batchFactor=0.8
hMax=( 0.03704 0.01235 0.00412 0.00138 0.00046 ) # 1/3^l ceiled with significance 1e-5

kernels=gengodunov # gengodunov or genmusclhancock

# Derived options

for fuseAlgorithmicSteps in "on" "off"
do 
  for patchSize in 7 11 15 19 # corresponds to orders=3 5 7 9
  do
    # SIMULATION END TIME
    T=( 0.01 0.00334 0.00112 0.00038 0.00013 )            # p=3
    if (( patchSize == 11 )); then
      T=( 0.006364 0.002126 0.000713 0.000242 0.000083 )  # p=5; (2*3+1)/(2*order+1)*T_3 ceiled with sig. 1e-6
    fi
    if (( patchSize == 15 )); then
      T=( 0.004667 0.001559 0.000523 0.000178 0.000061 )  # p=7
    fi
    if (( patchSize == 19 )); then
      T=( 0.003685 0.001231 0.000413 0.00014 0.000048 )   # p=9
    fi
    t=${T[i]}
    
    # Create script
    prefix=$project-$kernels
    if [ "$fuseAlgorithmicSteps" == "on" ]; then
      prefix+="-fused"
    else
      prefix+="-nonfused"
    fi
    script=convergence/hamilton.slurm-script
    newScript=convergence/hamilton-$prefix-N${patchSize}.slurm-script
    cp $script $newScript
   
    sed -i 's,prefix='$project',prefix='$prefix',g' $newScript
    sed -i 's,kernels=gen,kernels='$kernels',g' $newScript
    sed -i 's,p3,p'$patchSize',g' $newScript
    sed -i 's,regular-0,'$mesh',g' $newScript
    sed -i 's,script=convergence/hamilton.slurm-script,script='$newScript',g' $newScript
  
    for i in 0 1  2
    do
      mesh=regular-$i
      h=${hMax[i]}
      t=${T[i]}
       
      specPrefix=$prefix-$mesh
    
      # Create spec files
      coresPerTask=24
      spec=convergence/$project.exahype
      filename=convergence/$specPrefix-N$patchSize
      newSpec=$filename'.exahype'

      cp $spec $newSpec
      
      if [[ "$kernels" -eq "gengodunov" ]]; then
        sed -i -r 's,kernel(\s*)const(\s*)=(\s*)(.+),kernel\1const\2=\3generic::godunov,' $newSpec
      elif [[ "$kernels" -eq "genmuscl" ]]; then
        sed -i -r 's,kernel(\s*)const(\s*)=(\s*)(.+),kernel\1const\2=\3generic::musclhancock,' $newSpec
      fi

      sed -i -r 's,end-time(\s*)=(\s*)(([0-9]|\.)*),end-time\1=\2'$t',' $newSpec
      sed -i -r 's,ranks_per_node:([0-9]+),ranks_per_node:1,' $newSpec 
      sed -i -r 's,cores(\s+)=(\s+)([0-9]+),cores\1=\2'$coresPerTask',' $newSpec
     
      sed -i -r 's,skip-reduction-in-batched-time-steps(\s*)=(\s*)(\w+),skip-reduction-in-batched-time-steps\1=\2'$skipReductionInBatchedTimeSteps',' $newSpec
      sed -i -r 's,timestep-batch-factor(\s*)=(\s*)(([0-9]|\.)+),timestep-batch-factor\1=\2'$batchFactor',' $newSpec
      sed -i -r 's,fuse-algorithmic-steps(\s*)=(\s*)(\w+),fuse-algorithmic-steps\1=\2'$fuseAlgorithmicSteps',' $newSpec
    
      sed -i -r 's,patch-size(\s+)const(\s+)=(\s+)([0-9]+),patch-size\1const\2=\3'$patchSize',' $newSpec
      sed -i -r 's,maximum-mesh-size(\s*)=(\s*)(([0-9]|\.|-|\+|e|E)*),maximum-mesh-size\1=\2'$h',' $newSpec
    done
  done
done
