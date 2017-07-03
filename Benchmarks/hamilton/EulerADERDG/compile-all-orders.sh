exe1=ExaHyPE-EulerADERDG
exe2=ExaHyPE-EulerADERDG
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

  for p in 1 3 5 7 9
  do
    rm *.o
    sed -i -r "s,Order(\s+)= ([0-9]),Order\1= ${p}," $header
    make -j28 && \
    mv $exe1 $exe2-p$p-$SHAREDMEM-$COMPILER
#    sleep 2m
  done
done

