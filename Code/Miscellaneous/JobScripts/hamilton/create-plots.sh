python ../../../PerformanceAnalysis/plotRuntime.py -table x16-0.results.table x16-1.results.table x16-2.results.table x16-3.results.table x16-4.results.table x16-5.results.table -adapter ADERDGTimeStep Total -dimension 2 -xaxislabel "MPI ranks" -experimentdescription "h=0.7" "h=0.3" "h=0.1" "h=0.05" "h=0.01" "h=0.005" -output 16-ranks-per-node -singlecore 16

python ../../../PerformanceAnalysis/plotRuntime.py -table x8-0.results.table x8-1.results.table x8-2.results.table x8-3.results.table x8-4.results.table x8-5.results.table       -adapter ADERDGTimeStep Total -dimension 2 -xaxislabel "MPI ranks" -experimentdescription "h=0.7" "h=0.3" "h=0.1" "h=0.05" "h=0.01" "h=0.005" -output 8-ranks-per-node -singlecore 8

python ../../../PerformanceAnalysis/plotRuntime.py -table x4-0.results.table x4-1.results.table x4-2.results.table x4-3.results.table x4-4.results.table x4-5.results.table       -adapter ADERDGTimeStep Total -dimension 2 -xaxislabel "MPI ranks" -experimentdescription "h=0.7" "h=0.3" "h=0.1" "h=0.05" "h=0.01" "h=0.005" -output 4-ranks-per-node -singlecore 4



python ../../../PerformanceAnalysis/plotRuntime.py -table x16-0-noskip.results.table x16-1-noskip.results.table x16-2-noskip.results.table x16-3-noskip.results.table x16-4-noskip.results.table x16-5-noskip.results.table -adapter ADERDGTimeStep Total -dimension 2 -xaxislabel "MPI ranks" -experimentdescription "h=0.7" "h=0.3" "h=0.1" "h=0.05" "h=0.01" "h=0.005" -output 16-ranks-per-node-noskip -singlecore 16

python ../../../PerformanceAnalysis/plotRuntime.py -table x8-0-noskip.results.table x8-1-noskip.results.table x8-2-noskip.results.table x8-3-noskip.results.table x8-4-noskip.results.table x8-5-noskip.results.table       -adapter ADERDGTimeStep Total -dimension 2 -xaxislabel "MPI ranks" -experimentdescription "h=0.7" "h=0.3" "h=0.1" "h=0.05" "h=0.01" "h=0.005" -output 8-ranks-per-node-noskip -singlecore 8

python ../../../PerformanceAnalysis/plotRuntime.py -table x4-0-noskip.results.table x4-1-noskip.results.table x4-2-noskip.results.table x4-3-noskip.results.table x4-4-noskip.results.table x4-5-noskip.results.table       -adapter ADERDGTimeStep Total -dimension 2 -xaxislabel "MPI ranks" -experimentdescription "h=0.7" "h=0.3" "h=0.1" "h=0.05" "h=0.01" "h=0.005" -output 4-ranks-per-node-noskip -singlecore 4



#python ../../../PerformanceAnalysis/plotRuntime.py -table x16-0.results.table x16-1.results.table x16-2.results.table x16-3.results.table x16-4.results.table x16-5.results.table -adapter ADERDGTimeStep PredictorRerun Total -dimension 2 -xaxislabel "MPI ranks" -experimentdescription "h=0.7" "h=0.3" "h=0.1" "h=0.05" "h=0.01" "h=0.005" -output 16-ranks-per-node -singlecore 16

#python ../../../PerformanceAnalysis/plotRuntime.py -table x8-0.results.table x8-1.results.table x8-2.results.table x8-3.results.table x8-4.results.table x8-5.results.table -adapter ADERDGTimeStep PredictorRerun Total -dimension 2 -xaxislabel "MPI ranks" -experimentdescription "h=0.7" "h=0.3" "h=0.1" "h=0.05" "h=0.01" "h=0.005" -output 8-ranks-per-node -singlecore 8

#python ../../../PerformanceAnalysis/plotRuntime.py -table x4-0.results.table x4-1.results.table x4-2.results.table x4-3.results.table x4-4.results.table x4-5.results.table -adapter ADERDGTimeStep PredictorRerun Total -dimension 2 -xaxislabel "MPI ranks" -experimentdescription "h=0.7" "h=0.3" "h=0.1" "h=0.05" "h=0.01" "h=0.005" -output 4-ranks-per-node -singlecore 4

