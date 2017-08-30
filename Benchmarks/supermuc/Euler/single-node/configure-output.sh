echo "Configure project for single-node scaling test (output)."
mkdir single-node/results
rm -r *.o cfiles.mk ffiles.mk kernels
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/supermuc/Euler/single-node/Euler-output.exahype )
