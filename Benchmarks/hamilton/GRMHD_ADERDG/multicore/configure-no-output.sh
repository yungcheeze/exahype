echo "Configure project for multicore scaling test (no output)."
mkdir multicore/results
rm -rf *.o cfiles.mk ffiles.mk kernels cipofiles.mk
./projectpaths.cfg
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/GRMHD_ADERDG/multicore/GRMHD_ADERDG-no-output.exahype )

printf "\n\nPlease make sure you run the following command prior running make:\n"
echo "export COMPILER_LFLAGS=-lgsl -lgslcblas -lm"

