#!/bin/bash
echo "Building project ${PROJECT}..."
rm -f files.mk
make -f ${PROJECT}/makefile clean VERBOSE=1
