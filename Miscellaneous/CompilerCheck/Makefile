# This works on a typical 2017's Ubuntu installation with
# GCC 5.4, libtbb-dev 4.4, libopenmpi-dev 1.10

# It should also run without troubles with intel, MVAPICH
# and all that.

all: main-chrono main-tbb main-mpi main-all

# you want here g++ or icpc
CCX=icpc

main-chrono:
	$(CCX) -o $@ main-chrono.cc -std=c++11
	./$@
	rm main-chrono

main-tbb:
	$(CCX) -o $@ main-tbb.cc -ltbb
	./$@
	rm main-tbb

main-mpi:
	mpiCC -o $@ main-mpi.cc
	mpirun -np 4 ./$@
	rm main-mpi

main-all:
	mpiCC -o $@ main-all.cc -ltbb -std=c++11
	mpirun -np 4 ./$@
	rm main-all
