module load java

export EXAHYPE_CC="mpicc"
export COMPILER_CFLAGS="-DnoPackedRecords"

export MODE=Release
export COMPILER=Intel
export DISTRIBUTEDMEM=MPI
export ARCHITECTURE=hsw
export GPROF=off