source vars.sh

for m in 1 2
do
  if [ "$m" -eq "1" ]; then
    make clean
    export SHAREDMEM=TBB
  else
    make clean
    export SHAREDMEM=None
  fi

  echo "SHAREDMEM=$SHAREDMEM"
  read -p "press any key..."

  for p in 3 4 5 6 7 8 9
  do
    sed "s,int N=.*,int N=${p};," MyEulerSolver.cpp > MyEulerSolver.cpp_$p && \
    mv MyEulerSolver.cpp_${p} MyEulerSolver.cpp && \
    make -j4 && \
    mv ExaHyPE-Euler2d ExaHyPE-Euler2d-$SHAREDMEM-$p
	
#    sleep 2m
  done
done
