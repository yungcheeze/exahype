echo "Configure project for single-node scaling test (output)."
mkdir single-node/results
rm -rf *.o cfiles.mk ffiles.mk kernels cipofiles.mk
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/supermuc/Euler_FV/single-node/Euler_FV-output.exahype )
