## Preparation

It is important that we remove

PROJECT_CFLAGS+=-DnoMultipleThreadsMayTriggerMPICalls

from the local Makefile before we continue.




## Compile Code on Hamilton (Durham supercomputer)

module load intel/xe_2015.2
module load gcc/4.9.1

setenv SHAREDMEM none
setenv DISTRIBUTEDMEM none

make clean
make -j16
mv ExaHyPE-Z4 ExaHyPE-Z4-notbb-nompi




module load intel/xe_2015.2
module load gcc/4.9.1

setenv TBB_SHLIB "-L/gpfs/hamilton6/apps/intel/xe_2015/composer_xe_2015.2.164/tbb/lib/intel64/gcc4.4 -ltbb"

setenv SHAREDMEM tbb
setenv DISTRIBUTEDMEM none

make clean
make -j16
mv ExaHyPE-Z4 ExaHyPE-Z4-tbb-nompi




module load intel/xe_2015.2
module load gcc/4.9.1
module load intel_mpi

setenv SHAREDMEM none
setenv DISTRIBUTEDMEM mpi

make clean
make -j16
mv ExaHyPE-Z4 ExaHyPE-Z4-notbb-mpi




module load intel/xe_2015.2
module load gcc/4.9.1
module load intel_mpi
setenv TBB_SHLIB "-L/gpfs/hamilton6/apps/intel/xe_2015/composer_xe_2015.2.164/tbb/lib/intel64/gcc4.4 -ltbb"

setenv SHAREDMEM tbb
setenv DISTRIBUTEDMEM mpi
setenv PROJECT_CFLAGS

make clean
make -j16
mv ExaHyPE-Z4 ExaHyPE-Z4-tbb-mtmpi




module load intel/xe_2015.2
module load gcc/4.9.1
module load intel_mpi
setenv TBB_SHLIB "-L/gpfs/hamilton6/apps/intel/xe_2015/composer_xe_2015.2.164/tbb/lib/intel64/gcc4.4 -ltbb"

setenv SHAREDMEM tbb
setenv DISTRIBUTEDMEM mpi
setenv PROJECT_CFLAGS -DnoMultipleThreadsMayTriggerMPICalls

make clean
make -j16
mv ExaHyPE-Z4 ExaHyPE-Z4-tbb-mpi



Create link to experiments
ln -s ../../Benchmark_Mar2017 benchmarks
