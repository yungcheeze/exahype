module load java

export EXAHYPE_CC="mpicc"
export COMPILER_CFLAGS="-DnoParallelExchangePackedRecordsAtBoundary -DnoParallelExchangePackedRecordsBetweenMasterAndWorker -DnoParallelExchangePackedRecordsInHeaps -DnoParallelExchangePackedRecordsThroughoutJoinsAndForks"

export MODE=Release
export COMPILER=Intel
export DISTRIBUTEDMEM=MPI
export ARCHITECTURE=hsw
export GPROF=off

# optimised kernels
export USE_IPO=on

