echo "Configure project for multicore scaling test (output)."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/Euler_ADERDG/multicore/Euler_ADERDG-output.exahype )
mkdir multicore/results
rm *.o
