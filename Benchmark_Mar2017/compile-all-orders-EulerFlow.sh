exe1=ExaHyPE-Euler
exe2=ExaHyPE-EulerFlow
header=AbstractMyEulerSolver.h

sharedmemmode=notbb
for m in 1 2
do
  if [ "$m" -eq "1" ]; then
    make clean
    export SHAREDMEM=TBB
    sharedmemmode=tbb
  else
    make clean
    export SHAREDMEM=None
    sharedmemmode=notbb
  fi

  echo "SHAREDMEM=$SHAREDMEM"
  #read -p "press any key..."

  for p in 3 4 5 6 7 8 9
  do
    rm *.o
    sed -i -r "s,Order(\s+)= ([0-9]),Order\1= ${p}," $header
    make -j16 && \
    mv $exe1 $exe2'-p'$p'-'$sharedmemmode'-mpi'
#    sleep 2m
  done
done

