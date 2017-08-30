echo "Configure project for convergence test."
mkdir convergence/results
rm -r *.o cfiles.mk ffiles.mk kernels
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/Euler_FV/convergence/Euler_FV.exahype )
