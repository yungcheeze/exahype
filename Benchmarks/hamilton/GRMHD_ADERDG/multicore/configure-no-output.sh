echo "Configure project for multicore scaling test (no output)."
./projectpaths.cfg
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/GRMHD_ADERDG/multicore/GRMHD_ADERDG-no-output.exahype )
mkdir multicore/results
rm *.o

printf "\n\nPlease add the following line to your Makefile:\n"
echo "PROJECT_LFLAGS+=-lgsl -lgslcblas -lm"

