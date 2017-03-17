## Experiment description ##



 
Name                   | h    | Ranks per node    | threads  | Algorithmic variant
--------------------------------------------------------------------
EulerFlow-regular-0    | 0.01 | 28                | no tbb   | fused 
EulerFlow-regular-1    | 0.01 | 14                | 2        | fused
EulerFlow-regular-2    | 0.01 | 7                 | 4        | fused
EulerFlow-regular-3    | 0.01 | 4                 | 7        | fused
EulerFlow-regular-4    | 0.01 | 2                 | 14       | fused


All p=9 experiments are by 1/3 smaller than the p=3.

## Preparations

Compile executables and create spec files

```
cd ../ApplicationExamples/EulerFlow
ln -s ../../Benchmark_Mar2017 benchmarks

cp benchmarks/EulerFlow*.exahype .

source benchmarks/env_supermuc.sh
benchmarks/compile-all-orders-EulerFlow.sh
benchmarks/create-specfiles-for-all-orders.sh
```
Run the experiments

```
cp benchmarks/*.load-leveler # or *.slurm .

llsubmit EulerFlow... # or sbatch ...
```

### 

## Postprocessing of some output ##

python ~/workspace/peano/src/peano/performanceanalysis/domaindecompositionanalysis.py -dimension 3 -domainoffset -1 -1 -1 -domainsize 2 2 2 -file EulerFlow-28-ranks-per-node-8-nodes-0-p3-no-output-notbb.out

python ~/workspace/peano/src/peano/performanceanalysis/domaindecompositionanalysis.py -dimension 3 -domainoffset -1 -1 -1 -domainsize 2 2 2 -file EulerFlow-28-ranks-per-node-8-nodes-0-p9-no-output-notbb.out

