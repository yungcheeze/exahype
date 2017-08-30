echo "Configure project for convergence test."
mkdir convergence/results
rm -r *.o cfiles.mk ffiles.mk kernels
./projectpaths.cfg
( cd ../../../ && java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Benchmarks/hamilton/GRMHD_ADERDG/convergence/GRMHD_ADERDG.exahype )

printf "\n\nPlease make sure you run the following command prior running make:\n"
echo "export COMPILER_LFLAGS=-lgsl -lgslcblas -lm"
