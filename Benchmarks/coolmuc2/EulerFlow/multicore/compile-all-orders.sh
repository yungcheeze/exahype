
#path test case <-> root of Exahype
appcasepath=ApplicationExamples/EulerFlow
toroot=../..
#output or no-output
io=no-output

mpi=None

exe1=ExaHyPE-Euler
exe2=ExaHyPE-EulerFlow
compiler=Intel

for p in 9
do
  cd $toroot
  java -jar Toolkit/dist/ExaHyPE.jar --not-interactive $appcasepath/benchmarks/EulerFlow-$io-p$p-regular-0-t1-c1.exahype
  cd $appcasepath
  export DISTRIBUTEDMEM=$mpi
  
  for m in 1 2
  do
    if (( m == 1 )); then
      export SHAREDMEM=TBB
    else
      export SHAREDMEM=None
    fi

    echo "SHAREDMEM=$SHAREDMEM"
    #read -p "press any key..."
    
    make clean
    make -j56 && \
    mv $exe1 $exe2-p$p-$SHAREDMEM-$COMPILER
#    sleep 2m
  done
done

