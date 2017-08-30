echo "Configure project for multicore scaling test (no output)."
mkdir multicore/results
rm -r *.o cfiles.mk ffiles.mk kernels
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/Euler/multicore/Euler-no-output.exahype )