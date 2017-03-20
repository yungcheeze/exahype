#!/bin/bash

tasks_i=1

for i in 1 2 10 28
do 
  file=supermuc.load-leveler
  newfile=supermuc-$tasks_i-ranks-per-node-$i-nodes.load-leveler
  cp $file $newfile

  let tasks=tasks_i*i

  sed -i -r 's,NTASKS=([0-9]+),NTASKS='$tasks',g' $newfile
  sed -i -r 's,@ tasks_per_node(\s+)=(\s+)([0-9]+),@ tasks_per_node\1=\2'${tasks_i}',g' $newfile

  sed -i -r 's,@ node(\s+)=(\s+)([0-9]+),@ node\1=\2'$i',g' $newfile
  sed -i -r 's,SLURM_JOB_NUM_NODES=([0-9]+),SLURM_JOB_NUM_NODES='$i',g' $newfile

  sed -i 's,p3,p3,g' $newfile
done


#
sed -i -r 's,ranks_per_node:([0-9]+),ranks_per_node:'$tasks_i',g' EulerFlow-output.exahype
sed -i -r 's,ranks_per_node:([0-9]+),ranks_per_node:'$tasks_i',g' EulerFlow-no-output.exahype

./create-spec-for-all-orders.sh
