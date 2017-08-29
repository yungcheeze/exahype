echo "Configure project for single-node scaling test (output)."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/supermuc/Euler_ADERDG/single-node/Euler_ADERDG-output.exahype )
mkdir single-node/results
rm -r *.o cfiles.mk ffiles.mk kernels
