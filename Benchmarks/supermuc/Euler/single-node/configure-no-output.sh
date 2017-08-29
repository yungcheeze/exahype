echo "Configure project for single-node scaling test (no output)."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/supermuc/Euler/single-node/Euler-no-output.exahype )
mkdir single-node/results
rm -r *.o cfiles.mk ffiles.mk kernels
