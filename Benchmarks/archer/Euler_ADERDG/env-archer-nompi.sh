module load java/jdk1.8.0_51

module swap PrgEnv-cray PrgEnv-intel

export EXAHYPE_CC=icpc

export MODE=Release
export COMPILER=Intel
export DISTRIBUTEDMEM=None
export ARCHITECTURE=knl
export GPROF=off

# optimised kernels
export USE_IPO=on
export TBB_SHLIB="-L/opt/intel/compilers_and_libraries_2017.0.098/linux/tbb/lib/mic -ltbb"


