echo "Configure project for plenty-nodes scaling test (output)."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/supermuc/Euler/plenty-nodes/Euler-output.exahype )
mkdir plenty-nodes/results
rm *.o
