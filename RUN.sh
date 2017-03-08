#!/bin/sh

# stop on error
set -e

# download the grid generation framework Peano
./Peano/checkout-update-peano.sh

# compile the Java toolkit
./Toolkit/build.sh

# Optional: Install libxsmm

# git clone https://github.com/hfp/libxsmm.git Libxsmm
# cd Libxsmm/
# make generator

# Now you are ready to follow compile and run an ExaHyPE application
# according to the guidebook:
java -jar Toolkit/dist/ExaHyPE.jar --not-interactive ApplicationExamples/eulerflow2d.exahype

# set build parameters
export CC=Gnu
#export SHAREDMEM=TBB
#export TBB_INC=/usr/include/tbb
#export TBB_SHLIB="-L/usr/lib -ltbb"

# build sample application
cd ApplicationExamples/EulerFlow && make clean && make -j $(nproc)

# run sample application
./ExaHyPE-EulerFlow ../EulerFlow.exahype



