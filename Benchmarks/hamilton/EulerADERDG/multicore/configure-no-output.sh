echo "Configure project for multicore scaling test (no output)."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/EulerADERDG/multicore/EulerADERDG-no-output.exahype )
source env-hamilton-nompi.sh
