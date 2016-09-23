#!/bin/bash

MIN_DEGREE=3
MAX_DEGREE=8
NVAR=9

DEGREE=$MIN_DEGREE

cd src
while [ $DEGREE -le $MAX_DEGREE ]; do
	make clean
	make NVAR=$NVAR ORDER=$DEGREE EXE=main$DEGREE
	mv main$DEGREE ../
  let DEGREE=DEGREE+1;
done
cd ..
