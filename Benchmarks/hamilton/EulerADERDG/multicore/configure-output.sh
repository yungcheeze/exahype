echo "Configure project for multicore scaling test (output)."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/EulerADERDG/multicore/EulerADERDG-output.exahype )
mkdir multicore/results
