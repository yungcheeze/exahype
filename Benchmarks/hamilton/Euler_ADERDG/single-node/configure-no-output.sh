echo "Configure project for single-node scaling test (no output)."
mkdir single-node/results
rm -rf *.o cfiles.mk ffiles.mk kernels cipofiles.mk
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/Euler_ADERDG/single-node/Euler_ADERDG-no-output.exahype )
