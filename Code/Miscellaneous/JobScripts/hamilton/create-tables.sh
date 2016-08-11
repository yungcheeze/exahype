cp 1x1-0.results 1x4-0.results
cp 1x1-1.results 1x4-1.results
cp 1x1-2.results 1x4-2.results
cp 1x1-3.results 1x4-3.results
cp 1x1-4.results 1x4-4.results

cp 1x1-0.results 1x8-0.results
cp 1x1-1.results 1x8-1.results
cp 1x1-2.results 1x8-2.results
cp 1x1-3.results 1x8-3.results
cp 1x1-4.results 1x8-4.results

cp 1x1-0.results 1x16-0.results
cp 1x1-1.results 1x16-1.results
cp 1x1-2.results 1x16-2.results
cp 1x1-3.results 1x16-3.results
cp 1x1-4.results 1x16-4.results

cp 2x2-0.results 2x4-0.results
cp 2x2-1.results 2x4-1.results
cp 2x2-2.results 2x4-2.results
cp 2x2-3.results 2x4-3.results
cp 2x2-4.results 2x4-4.results

cp 2x2-0.results 2x8-0.results
cp 2x2-1.results 2x8-1.results
cp 2x2-2.results 2x8-2.results
cp 2x2-3.results 2x8-3.results
cp 2x2-4.results 2x8-4.results

cp 2x2-0.results 2x16-0.results
cp 2x2-1.results 2x16-1.results
cp 2x2-2.results 2x16-2.results
cp 2x2-3.results 2x16-3.results
cp 2x2-4.results 2x16-4.results

python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x4-0.results -maxnodes 1025
python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x4-1.results -maxnodes 1025
python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x4-2.results -maxnodes 1025
python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x4-3.results -maxnodes 1025
python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x4-4.results -maxnodes 1025
python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x4-5.results -maxnodes 1025


python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x8-0.results -maxnodes 1025
python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x8-1.results -maxnodes 1025
python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x8-2.results -maxnodes 1025
python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x8-3.results -maxnodes 1025
python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x8-4.results -maxnodes 1025
python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x8-5.results -maxnodes 1025


python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x16-0.results -maxnodes 1025
python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x16-1.results -maxnodes 1025
python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x16-2.results -maxnodes 1025
python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x16-3.results -maxnodes 1025
python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x16-4.results -maxnodes 1025
python ../../../PerformanceAnalysis/extractTable.py -path . -prefix "" -postfix x16-5.results -maxnodes 1025



