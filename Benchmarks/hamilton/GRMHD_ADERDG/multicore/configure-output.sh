echo "Configure project for multicore scaling test (output)."
mkdir multicore/results
rm -r *.o cfiles.mk ffiles.mk kernels
./projectpaths.cfg
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/GRMHD_ADERDG/multicore/GRMHD_ADERDG-output.exahype )

printf "\n\nPlease make sure you run the following command prior running make:\n"
echo "export COMPILER_LFLAGS=-lgsl -lgslcblas -lm"
