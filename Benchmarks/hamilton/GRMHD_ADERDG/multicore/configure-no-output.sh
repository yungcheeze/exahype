echo "Configure project for multicore scaling test (no output)."
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/GRMHD_ADERDG/multicore/GRMHD_ADERDG-no-output.exahype )
mkdir multicore/results
