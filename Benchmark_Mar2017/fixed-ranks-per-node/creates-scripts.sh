#!/bin/bash

tasks_i=1

for i in 1 2 4 8 16
do 
  file=supermuc.load-leveler
  newfile=supermuc-$tasks_i-ranks-per-node-$i-nodes.load-leveler
  cp $file $newfile

  let tasks=tasks_i*i

  sed -i -r 's,NTASKS=([0-9]+),NTASKS='$tasks',g' $newfile
  sed -i -r 's,@ tasks_per_node(\s+)=(\s+)([0-9]+),@ tasks_per_node\1=\2'${tasks_i}',g' $newfile

  sed -i -r 's,@ node(\s+)=(\s+)([0-9]+),@ node\1=\2'$i',g' $newfile

  sed -i 's,p3,p3,g' $newfile

done
