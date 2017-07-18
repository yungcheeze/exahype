exe1=ExaHyPE-Euler_FV
exe2=ExaHyPE-Euler_FV
header=AbstractEulerSolver.h


for m in 1 2
do
  if (( m == 1 )); then
    make clean
    export SHAREDMEM=TBB
  else
    make clean
    export SHAREDMEM=None
  fi

  echo "SHAREDMEM=$SHAREDMEM"
  #read -p "press any key..."

  for N in 7 11 15 17 # corresponds to orders 3 5 7 9
  do
    rm *.o
    sed -i -r "s,PatchSize(\s+)= ([0-9]+),PatchSize\1= ${N}," $header
    make -j28 && \
    mv $exe1 $exe2-p$p-$SHAREDMEM-$COMPILER
#    sleep 2m
  done
done

