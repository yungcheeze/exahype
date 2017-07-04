echo "Configure project for convergence test."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/EulerFV/convergence/EulerFV.exahype )
mkdir convergence/results
