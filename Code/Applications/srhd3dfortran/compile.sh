#!/bin/sh

set -e

cd ../../
java -jar ExaHyPE.jar  --not-interactive Applications/srhd3dfortran.exahype

cd -
make clean
make -j $(nproc)
