# unpack the grid generation framework Peano
cd Code/Peano
tar xvfz peano.tar.gz
git checkout .gitignore
cd ..

# run toolkit for the generation of user applications
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


