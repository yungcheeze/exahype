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
  #read -p "press any key..."

  for p in 3 4 5 6 7 8 9
  do
    rm *.o

    let N=2*p+1
    echo $N
    sed -i -E "s,Order(\s+)= ([0-9]+),Order\1= ${p}," AbstractLimitingADERDG_ADERDG.h
    sed -i -E "s,PatchSize(\s+)= ([0-9]+),PatchSize\1= ${N}," AbstractLimitingADERDG_FV.h
    make -j4 && \
    mv ExaHyPE-Euler ExaHyPE-Euler-$SHAREDMEM-$p

#    sleep 2m
  done
done

