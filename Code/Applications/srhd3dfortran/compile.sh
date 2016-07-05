#!/bin/sh

export CC=gcc
export SHAREDMEM=TBB
export TBB_INC=/usr/include/tbb
export TBB_SHLIB="-L/usr/lib -ltbb"

# enables -g -O0 etc.
#export MODE="Debug"


set -e

cd ../../
java -jar ExaHyPE.jar  --not-interactive Applications/srhd3dfortran.exahype

cd -
#make clean
make -j $(nproc)
