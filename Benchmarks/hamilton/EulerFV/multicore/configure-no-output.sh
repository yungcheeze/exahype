echo "Configure project for multicore scaling test (no output)."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/EulerFV/multicore/EulerFV-no-output.exahype )
mkdir multicore/results
