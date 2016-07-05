#!/bin/bash

cd $(dirname "$0")

# Vasco, 26. June by mail: No Sharedmem in this user-spec file.
export CC=gcc
export SHAREDMEM=None
#export TBB_INC=/usr/include/tbb

MPI_LDFLAGS="$(mpicc -showme:link)"
# at ubuntu 14: -pthread -L/usr//lib -L/usr/lib/openmpi/lib -lmpi -ldl -lhwloc

export TBB_SHLIB="-L/usr/lib -ltbb $MPI_LDFLAGS"

set -e

cd ../../

# Vasco: do not run toolkit on SRHD right now, as this is the new Fortran prototype.

#java -jar ExaHyPE.jar  --not-interactive Applications/srhd3dfortran.exahype

cd -
#make clean
make -j $(nproc)

