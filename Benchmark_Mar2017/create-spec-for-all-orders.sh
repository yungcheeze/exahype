
for specfile in 'EulerFlow-output' 'EulerFlow-no-output'
do
  for p in 3 4 5 6 7 8 9
  do
    cp $specfile'.exahype' $specfile'_'$p.'exahype'
    sed -i -r "s,order(\s*)const(\s*)=(\s*)([0-9]),order\1const\2=\3${p}," $specfile'_'$p.'exahype'
  done
done

