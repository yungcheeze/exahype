MAIN=mainC

echo "Clean"
rm -f *.o
rm -f $MAIN

echo "Compile"
icpc -c -std=c++11 -O3 -DALIGNMENT=32 -xHost -restrict DGMatrices.cpp
icpc -c -std=c++11  -O3 -DALIGNMENT=32 -xHost -restrict GaussLegendreQuadrature.cpp
icpc -c -std=c++11  -O3 -DALIGNMENT=32 -xHost -restrict volumeIntegral.cpp
icpc -c -std=c++11  -O3 -DALIGNMENT=32 -xHost -restrict main.cpp
icpc -std=c++11  -O3 -DALIGNMENT=32 -xHost -restrict DGMatrices.o GaussLegendreQuadrature.o volumeIntegral.o main.o -o $MAIN
echo "Done"
