#!/bin/bash

function verbose { echo $@; $@; }

set -e

verbose g++ -c -g3 --std=c++11 -DDim3 -Wall -I../../Peano/ \
	indi-test.cpp InitialData.cpp \
	PDE/Cons2Prim.cpp PDE/Prim2Cons.cpp PDE/PDE.cpp 
	
verbose g++ -g3 -oinditest.out \
	FortranInitialData/InitialDataFort.o Fortran/*.o \
	*.o -lgfortran
