echo "Configure project for multicore scaling test (no output)."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/Euler_ADERDG/multicore/Euler_ADERDG-no-output.exahype )
mkdir multicore/results
rm -r *.o cfiles.mk ffiles.mk kernels
