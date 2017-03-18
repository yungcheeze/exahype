#!/bin/bash
for i in 1 2 4 7 8 14 16 28
do
  sed -r 's,NTASKS=([0-9]+),NTASKS='$i',g' supermuc-1-nodes.load-leveler > supermuc-$i-ranks-per-node-1-nodes.load-leveler
done
