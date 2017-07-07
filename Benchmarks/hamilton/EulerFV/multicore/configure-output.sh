echo "Configure project for multicore scaling test (output)."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/EulerFV/multicore/EulerFV-output.exahype )
mkdir multicore/results
