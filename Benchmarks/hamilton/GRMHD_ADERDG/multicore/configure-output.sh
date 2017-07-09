echo "Configure project for multicore scaling test (output)."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/GRMHD_ADERDG/multicore/GRMHD_ADERDG-output.exahype )
mkdir multicore/results
