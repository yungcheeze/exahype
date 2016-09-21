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
CFLAGS=-I. -std=c++11 -Wall -restrict -O2 -xHost
PREPVAR = -DNDEBUG -DALIGNMENT=64 -DFILESUFFIX=O2 -DNVAR=9 -DNVARPAD=16
MAIN=test_O2



OBJ = 	./glue/DGMatricesGeneric.o \
		./srcGen/DGMatrices3.o \
		./srcGen/DGMatrices4.o \
		./srcGen/DGMatrices5.o \
		./srcGen/DGMatrices6.o \
		./srcGen/DGMatrices7.o \
		./srcGen/DGMatrices8.o \
		./srcGen/GaussLegendreQuadrature3.o \
		./srcGen/GaussLegendreQuadrature4.o \
		./srcGen/GaussLegendreQuadrature5.o \
		./srcGen/GaussLegendreQuadrature6.o \
		./srcGen/GaussLegendreQuadrature7.o \
		./srcGen/GaussLegendreQuadrature8.o \
		./srcGen/volumeIntegral3.o \
		./srcGen/volumeIntegral4.o \
		./srcGen/volumeIntegral5.o \
		./srcGen/volumeIntegral6.o \
		./srcGen/volumeIntegral7.o \
		./srcGen/volumeIntegral8.o \
		./glue/volumeKernelTest.o

.PHONY: default
default: all

%.o : %.cpp
	$(CC) -c $(CFLAGS) $(PREPVAR) $*.cpp -o $*.o

all: $(OBJ)
	$(CC) $(CFLAGS) $(PREPVAR) -o $(MAIN) $(OBJ)
	
clean:
	rm $(OBJ) $(MAIN)
