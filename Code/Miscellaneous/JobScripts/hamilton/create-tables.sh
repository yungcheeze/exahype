for experiment in 0 1 2 3 4 5
do
  for extension in "" "-batch-noskip-exchange" "-batch-noskip-noexchange" "-batch-skip-exchange" "-batch-skip-noexchange" "-noexchange"
  do
    cp 1x1-$experiment.results 1x4-$experiment$extension.results
    cp 1x1-$experiment.results 1x8-$experiment$extension.results
    cp 1x1-$experiment.results 1x16-$experiment$extension.results
  done
done


for experiment in 0 1 2 3 4 5
do
  for extension in "" "-batch-noskip-exchange" "-batch-noskip-noexchange" "-batch-skip-exchange" "-batch-skip-noexchange" "-noexchange"
  do
    python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x4-$experiment$extension.results -maxnodes 5000
    python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x8-$experiment$extension.results -maxnodes 5000
    python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x16-$experiment$extension.results -maxnodes 5000
  done
done



