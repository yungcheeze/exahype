echo "Configure project for convergence test."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/supermuc/Euler_ADERDG/convergence/Euler_ADERDG.exahype )
mkdir convergence/results
rm -r *.o cfiles.mk ffiles.mk kernels
