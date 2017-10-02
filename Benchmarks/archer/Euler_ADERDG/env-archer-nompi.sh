module load java/jdk1.8.0_51

module swap PrgEnv-cray PrgEnv-intel
module load gcc

export EXAHYPE_CC=CC

export MODE=Release
export COMPILER=Intel
#export COMPILER=gnu
#export COMPILER=Cray
#export DISTRIBUTEDMEM=None
export ARCHITECTURE=knl
export GPROF=off

# optimised kernels
export USE_IPO=on
#export USE_IPO=off
export TBB_SHLIB="-L/opt/intel/compilers_and_libraries_2017.0.098/linux/tbb/lib/intel64/gcc4.7 -ltbb"

# Tell the Cray environment that we wanna use dynamic libraries (for TBB)
CRAYPE_LINK_TYPE=dynamic


