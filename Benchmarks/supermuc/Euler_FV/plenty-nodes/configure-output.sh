echo "Configure project for plenty-nodes scaling test (output)."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/supermuc/Euler_FV/plenty-nodes/Euler_FV-output.exahype )
mkdir plenty-nodes/results
rm -r *.o cfiles.mk ffiles.mk kernels
