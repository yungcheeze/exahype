for experiment in 0 1 2 3 4 5
do
  cp 1x1-$experiment.results 1x4-$experiment.results
  cp 1x1-$experiment.results 1x8-$experiment.results
  cp 1x1-$experiment.results 1x16-$experiment.results

  cp 2x2-$experiment.results 2x4-$experiment.results
  cp 2x2-$experiment.results 2x8-$experiment.results
  cp 2x2-$experiment.results 2x16-$experiment.results

  cp 1x1-$experiment-noskip.results 1x4-$experiment-noskip.results
  cp 1x1-$experiment-noskip.results 1x8-$experiment-noskip.results
  cp 1x1-$experiment-noskip.results 1x16-$experiment-noskip.results

  cp 2x2-$experiment-noskip.results 2x4-$experiment-noskip.results
  cp 2x2-$experiment-noskip.results 2x8-$experiment-noskip.results
  cp 2x2-$experiment-noskip.results 2x16-$experiment-noskip.results
done


for experiment in 0 1 2 3 4 5
do
  python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x4-$experiment.results -maxnodes 1025
  python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x8-$experiment.results -maxnodes 1025
  python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x16-$experiment.results -maxnodes 1025

  python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x4-$experiment-noskip.results -maxnodes 1025
  python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x8-$experiment-noskip.results -maxnodes 1025
  python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x16-$experiment-noskip.results -maxnodes 1025
done



