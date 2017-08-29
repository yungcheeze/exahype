echo "Configure project for plenty-nodes scaling test (no output)."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/supermuc/Euler_ADERDG/plenty-nodes/Euler_ADERDG-no-output.exahype )
mkdir plenty-nodes/results
rm -r *.o cfiles.mk ffiles.mk kernels
