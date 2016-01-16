#!/bin/bash

echo Postprocess terminal output of the OracleForOnePhaseWithGrainSizeSampling
echo \(C\) 2013 Tobias Weinzierl
echo Usage: anydir/postprocess-grain-size-sampling.sh output-file adapter-number
echo "       where the script postprocess-grain-size-sampling.sh is in the same" 
echo "       directory as the Python script with the same name and its corresponding"
echo "       gnuplot instruction file "

script_dir=$(dirname $0)

echo Assume all scripts to reside in $script_dir

python $script_dir/postprocess-grain-size-sampling.py $1 $2 $script_dir
