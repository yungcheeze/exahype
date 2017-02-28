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
    rm *.o
    sed -i -E "s,Order(\s+)= ([0-9]),Order\1= ${p}," AbstractMyEulerSolver.h
    make -j4 && \
    mv ExaHyPE-Euler ExaHyPE-Euler-$SHAREDMEM-$p

#    sleep 2m
  done
done

