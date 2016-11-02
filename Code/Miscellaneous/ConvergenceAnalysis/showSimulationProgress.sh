#!/bin/bash
# shows a listing about the proress of all simulations
# Output is a valid CSV file.

h=$(mktemp)
SEP="\t"

echo -en "SimulationName${SEP}"
echo -en "NumReductionFiles${SEP}"
echo -en "EachRedLength${SEP}"
echo -en "NumLinesLog${SEP}"
echo -en "FinishedStatus${SEP}"
echo -e  "Walltime"

for sim in $(find simulations/ -maxdepth 1 -type d | sort); do
	[[ $sim == "simulations/" ]] && continue
	printf "%-40s" $sim
	if ls $sim/output/*.asc > $h 2>/dev/null; then
		echo -ne "$SEP" $(ls $(<$h) | wc -l) # number of files
		echo -ne "$SEP" $(wc -l $(<$h) | head -n-1 | awk '{ print $1 }' | uniq) # entries in each file
		echo -ne "$SEP" $(wc -l $sim/*log | head -n1 | awk '{ print $1 }') # number of lines in logfile
		if T=$(grep -oh "Finished ExaHyPE successfully" $sim/*.log); then
			echo -ne "$SEP FINISHED"
		elif T=$(grep -oh "ExaHyPE binary failed" $sim/*.log); then
			echo -ne "$SEP FAILED"
		else
			echo -ne "$SEP '" $(grep -h step $sim/*.log | tail -n1 | awk '{print $6" "$7" "$8" "$9 }') "'"
		fi
		#if T=$(tail $sim/*.log | grep user | awk '{ print $4 }'); then
		if T=$(grep -E "user\s+[0-9]" $sim/*.log | awk '{ print $4 }'); then
			echo -ne "$SEP $T"
		else	echo -ne "$SEP noWalltime"
		fi
	else
		echo -n "noResults"
	fi
	echo # newline 
done 
