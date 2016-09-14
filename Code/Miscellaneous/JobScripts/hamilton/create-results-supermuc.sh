#!/bin/bash

#rm *.results


MyNODES=1
MyNTASKS_PER_NODE=1
cp euler-test.load-leveler tmp.load-leveler
sed -i "s/NODES/$MyNODES/g" tmp.load-leveler
sed -i "s/NTASKS_PER_NODE/$MyNTASKS_PER_NODE/g" tmp.load-leveler
#llsubmit tmp.load-leveler


for MyNODES in 1 2 3 4 5 6 7 8 12 16 20 24 28 32
do
 for MyNTASKS_PER_NODE in 4 8 16
 do 
  cp euler-test.load-leveler tmp.load-leveler
  sed -i "s/NODES/$MyNODES/g" tmp.load-leveler
  sed -i "s/NTASKS_PER_NODE/$MyNTASKS_PER_NODE/g" tmp.load-leveler
  #llsubmit tmp.load-leveler
 done
done

for MyNODES in 40 48 56 64 80 96 112 128 
do
 for MyNTASKS_PER_NODE in 4 8 16
 do
  cp euler-general.load-leveler tmp.load-leveler
  sed -i "s/NODES/$MyNODES/g" tmp.load-leveler
  sed -i "s/NTASKS_PER_NODE/$MyNTASKS_PER_NODE/g" tmp.load-leveler
  #llsubmit tmp.load-leveler
 done
done


llq -u $USER

