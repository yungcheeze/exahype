source /etc/profile.d/modules.sh

module purge
module load java
module load slurm
module load intel/xe_2017.2
module load intelmpi/intel/2017.2
module load gcc

export EXAHYPE_CC=mpicc

export MODE=Release
export COMPILER=Intel
export DISTRIBUTEDMEM=MPI
