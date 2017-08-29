echo "Configure project for multicore scaling test (no output)."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/Euler/multicore/Euler-no-output.exahype )
mkdir multicore/results
rm -r *.o cfiles.mk ffiles.mk kernels