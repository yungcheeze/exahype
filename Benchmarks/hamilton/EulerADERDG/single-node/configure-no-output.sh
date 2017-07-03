echo "Configure project for single-node scaling test (no output)."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/EulerADERDG/single-node/EulerADERDG-no-output.exahype )
source ../env-hamilton-mpi.sh
