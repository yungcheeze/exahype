echo "Configure project for single-node scaling test (no output)."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/Euler_FV/single-node/Euler_FV-no-output.exahype )
mkdir single-node/results
rm -r *.o cfiles.mk ffiles.mk kernels
