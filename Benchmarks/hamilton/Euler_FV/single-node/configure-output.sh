echo "Configure project for single-node scaling test (output)."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/Euler_FV/single-node/Euler_FV-output.exahype )
mkdir single-node/results
rm -r *.o cfiles.mk ffiles.mk kernels
