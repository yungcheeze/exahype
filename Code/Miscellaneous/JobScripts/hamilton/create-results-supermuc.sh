#!/bin/bash

#rm *.results

for MyExperiment in 0 1 2 3 4
do
  MyNODES=1
  MyNTASKS_PER_NODE=1
  Script=tmp.$MyNODES.$MyNTASKS_PER_NODE.$MyExperiment.load-leveler
  cp euler.load-leveler $Script
  sed -i "s/NODES/$MyNODES/g" $Script
  sed -i "s/NTASKS_PER_NODE/$MyNTASKS_PER_NODE/g" $Script
  sed -i "s/EXPERIMENT/$MyExperiment/g" $Script
  llsubmit $Script
 
  for MyNODES in 1 2 3 4 5 6 7 8
  do
   for MyNTASKS_PER_NODE in 4 8 16
   do 
    Script=tmp.$MyNODES.$MyNTASKS_PER_NODE.$MyExperiment.load-leveler
    cp euler.load-leveler $Script
    sed -i "s/NODES/$MyNODES/g" $Script
    sed -i "s/NTASKS_PER_NODE/$MyNTASKS_PER_NODE/g" $Script
    sed -i "s/EXPERIMENT/$MyExperiment/g" $Script
    llsubmit $Script
   done
  done
done

llq -u $USER

