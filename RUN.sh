# unpack the grid generation framework Peano
cd Code/Peano
tar xvfz peano.tar.gz
git checkout .gitignore
cd ..

# run toolkit for the generation of user applications
java -jar Toolkit/dist/ExaHyPE.jar --not-interactive Applications/eulerflow2d.exahype

# set build parameters
export CC=gcc
#export SHAREDMEM=TBB
#export TBB_INC=/usr/include/tbb
#export TBB_SHLIB="-L/usr/lib -ltbb"

# build sample application
cd Applications/eulerflow2d && make -j

# run sample application
./ExaHyPE-Euler2d ../eulerflow2d.exahype


