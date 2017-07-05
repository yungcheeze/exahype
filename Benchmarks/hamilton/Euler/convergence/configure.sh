rm *.o
echo "Configure project for convergence test."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/Euler/convergence/Euler.exahype )
mkdir convergence/results
