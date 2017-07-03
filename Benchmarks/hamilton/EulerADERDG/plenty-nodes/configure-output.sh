echo "Configure project for plenty-nodes scaling test (output)."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/EulerADERDG/plenty-nodes/EulerADERDG-output.exahype )
source ../env-hamilton-mpi.sh
