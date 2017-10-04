echo "Configure project for multicore scaling test (no output)."
mkdir multicore/results
rm -r *.o cfiles.mk ffiles.mk kernels
( cd ../../../ && java -XX:MaxHeapSize=512m -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/archer/Euler_ADERDG/multicore/Euler_ADERDG-no-output.exahype )
