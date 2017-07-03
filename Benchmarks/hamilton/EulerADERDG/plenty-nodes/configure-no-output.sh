echo "Configure project for plenty-nodes scaling test (no output)."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/EulerADERDG/plenty-nodes/EulerADERDG-no-output.exahype )
source ../env-hamilton-mpi.sh
