#!/bin/sh

set -e

cd ../../
java -jar ExaHyPE.jar  --not-interactive Applications/srhd2dcpp.exahype

cd -
make clean
make -j $(nproc)
