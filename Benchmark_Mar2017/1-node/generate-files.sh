#!/bin/bash

for i in 1 2 10 28
do 
  #@ tasks_per_node = 28
  file=supermuc.load-leveler
  newfile=supermuc-$i-ranks-per-node-1-nodes.load-leveler
  cp $file $newfile

  sed -i -r 's,NTASKS=([0-9]+),NTASKS='$i',g' $newfile
  sed -i -r 's,tasks_per_node(\s+)=(\s+)([0-9]+),tasks_per_node\1=\2'$i',g' $newfile
  sed -i    's,p3,p3,g' $newfile

  #
  file=EulerFlow-output.exahype
  prefix1=EulerFlow-output-$i-ranks-per-node
  newfile=$prefix1'.exahype'
  cp $file $newfile

  sed -i -r 's,ranks_per_node:([0-9]+),ranks_per_node:'$i',g' $newfile

  
  file=EulerFlow-no-output.exahype
  prefix2=EulerFlow-no-output-$i-ranks-per-node
  newfile=$prefix2.exahype
  cp $file $newfile

  sed -i -r 's,ranks_per_node:([0-9]+),ranks_per_node:'$i',g' $newfile

  ./create-spec-for-all-orders.sh $prefix1 $prefix2
done
