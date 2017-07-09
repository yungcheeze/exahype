echo "Configure project for multicore scaling test (no output)."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/Euler_FV/multicore/Euler_FV-no-output.exahype )
mkdir multicore/results
rm *.o
