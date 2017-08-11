echo "Configure project for plenty-nodes scaling test (no output)."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/Euler_FV/plenty-nodes/Euler_FV-no-output.exahype )
mkdir plenty-nodes/results
rm *.o
