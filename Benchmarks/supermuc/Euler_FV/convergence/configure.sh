echo "Configure project for convergence test."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/supermuc/Euler_FV/convergence/Euler_FV.exahype )
mkdir convergence/results
rm *.o
