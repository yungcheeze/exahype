# This file is part of the ExaHyPE project.
# Copyright (c) 2016  http://exahype.eu
# All rights reserved.
# 
# The project has received funding from the European Union's Horizon 
# 2020 research and innovation programme under grant agreement
# No 671698. For copyrights and licensing, please consult the webpage.
# 
# Released under the BSD 3 Open Source License.
# For the full license text, see LICENSE.txt


CC=icpc
CFLAGS=-I. -std=c++11 -Wall -restrict -O2 -xCORE-AVX2
PREPVAR = -DNDEBUG -DALIGNMENT=64 -DFILESUFFIX=O2_9_hsw -DNVAR=9 -DNVARPAD=12
MAIN=test_O2_9_hsw



OBJ = 	./srcGen/GaussLegendreQuadrature3.o \
		./srcGen/GaussLegendreQuadrature4.o \
		./srcGen/GaussLegendreQuadrature5.o \
		./srcGen/GaussLegendreQuadrature6.o \
		./srcGen/GaussLegendreQuadrature7.o \
		./srcGen/GaussLegendreQuadrature8.o \
		./srcGen/predictor3.o \
		./srcGen/predictor4.o \
		./srcGen/predictor5.o \
		./srcGen/predictor6.o \
		./srcGen/predictor7.o \
		./srcGen/predictor8.o \
		./glue/predictorKernelTest.o

.PHONY: default
default: all

%.o : %.cpp
	$(CC) -c $(CFLAGS) $(PREPVAR) $*.cpp -o $*.o

all: $(OBJ)
	$(CC) $(CFLAGS) $(PREPVAR) -o $(MAIN) $(OBJ)
	
clean:
	rm $(OBJ) $(MAIN)
