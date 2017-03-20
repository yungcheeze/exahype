#!/bin/bash

sed -i -r 's,TASK_ID=([0-9]+),TASK_ID=1,g' supermuc-1-nodes.load-leveler

for i in 1 2 4 7 8 14 16 28
do 
  #@ tasks_per_node = 28i
  newfile=supermuc-$i-ranks-per-node-1-nodes.load-leveler

  sed -r 's,NTASKS=([0-9]+),NTASKS='$i',g' supermuc-1-nodes.load-leveler > $newfile
  sed -i -r 's,tasks_per_node(\s+)=(\s+)([0-9]+),tasks_per_node\1=\2'$i',g' $newfile
  sed -i 's,p3,p9,g' $newfile
done
