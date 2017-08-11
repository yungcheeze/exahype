echo "Configure project for plenty-nodes scaling test (no output)."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/Euler/plenty-nodes/Euler-no-output.exahype )
mkdir plenty-nodes/results
rm *.o
