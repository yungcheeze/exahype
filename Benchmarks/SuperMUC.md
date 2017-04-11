## Preparation

It is important that we remove

PROJECT_CFLAGS+=-DnoMultipleThreadsMayTriggerMPICalls

from the local Makefile	before we continue. This is a default flag that we switch on/off here manually.

For all the experiments, we use IBM MPI. This is the same variant the SeisSol 
team uses for their runs according to Sebastian Rettenberger. 



## Compile Code on SuperMUC

# Regenerate all the stuff
java -jar Toolkit/dist/ExaHyPE.jar --not-interactive AstroApplications/Z4/benchmarks/Z4-regular-0-p3-output.exahype
vi AstroApplications/Z4/Makefile
export ORDER=3


# or
java -jar Toolkit/dist/ExaHyPE.jar --not-interactive AstroApplications/Z4/benchmarks/Z4-regular-0-p9-output.exahype
vi AstroApplications/Z4/Makefile
export ORDER=9


# Trigger builds
module load gcc
export SHAREDMEM=none
export DISTRIBUTEDMEM=none

make clean
make -j16
mv ExaHyPE-Z4 ExaHyPE-Z4-p$ORDER-notbb-nompi




module load gcc
module load tbb
export SHAREDMEM=tbb
export DISTRIBUTEDMEM=none

make clean
make -j16
mv ExaHyPE-Z4 ExaHyPE-Z4-p$ORDER-tbb-nompi



module load gcc
module load tbb
export SHAREDMEM=none
export DISTRIBUTEDMEM=mpi
export EXAHYPE_CC=mpiCC
make clean
make -j16
mv ExaHyPE-Z4 ExaHyPE-Z4-p$ORDER-notbb-mpi




module load gcc
module load tbb
export SHAREDMEM=tbb
export DISTRIBUTEDMEM=mpi
export EXAHYPE_CC=mpiCC

make clean
make -j16
mv ExaHyPE-Z4 ExaHyPE-Z4-p$ORDER-tbb-mtmpi




module load gcc
module load tbb
export SHAREDMEM=tbb
export DISTRIBUTEDMEM=mpi
export EXAHYPE_CC=mpiCC
export PROJECT_CFLAGS=-DnoMultipleThreadsMayTriggerMPICalls

make clean
make -j16
mv ExaHyPE-Z4 ExaHyPE-Z4-p$ORDER-tbb-mpi



module load gcc
module load tbb
export SHAREDMEM=tbb
export DISTRIBUTEDMEM=mpi
export EXAHYPE_CC=mpiCC
export PROJECT_CFLAGS=-DnoMultipleThreadsMayTriggerMPICalls
export MODE=Profile

make clean
make -j16
mv ExaHyPE-Z4 ExaHyPE-Z4-p$ORDER-tbb-mpi-profile





module load gcc
module load tbb
module load scalasca
export SHAREDMEM=tbb
export DISTRIBUTEDMEM=mpi
export EXAHYPE_CC="scalasca -instrument mpiCC"
export EXAHYPE_FC="scalasca -instrument ifort"
export PROJECT_CFLAGS=-DnoMultipleThreadsMayTriggerMPICalls -g3
export MODE=Release

make clean
make -j16
mv ExaHyPE-Z4 ExaHyPE-Z4-p$ORDER-tbb-mpi-scalasca




## Prepare experiment environment

Create link to experiments
ln -s ../../Benchmark_Mar2017 benchmarks

rm *.properties

cp benchmarks/exahype.log-filter .

mkdir global-integrals
mkdir output
