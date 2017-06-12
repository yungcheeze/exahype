
#path test case <-> root of Exahype
appcasepath=ApplicationExamples/EulerFlow
toroot=../..
#output or no-output
io=no-output

exe1=ExaHyPE-Euler
exe2=ExaHyPE-EulerFlow


for m in 1 2
do
  for p in 9
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

    cd $toroot
    java -jar Toolkit/dist/ExaHyPE.jar $appcasepath/benchmarks/EulerFlow-$io-p$p-regular-0-t1-c1.exahype
    cd $appcasepath
    
    make -j56 && \
    mv $exe1 $exe2-p$p-$SHAREDMEM-$COMPILER
#    sleep 2m
  done
done

