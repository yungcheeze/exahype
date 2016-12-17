#!/bin/bash
# shows a listing about the proress of all simulations
# Output is a valid CSV file.

h=$(mktemp)
SEP="\t"

# where to look for simulations
SIMBASE="${SIMBASE:=simulations/}"

echo -en "SimulationName${SEP}"
echo -en "NumReductionFiles${SEP}"
echo -en "EachRedLength${SEP}"
echo -en "NumLinesLog${SEP}"
echo -en "FinishedStatus${SEP}"
echo -ne "LastTimeStep${SEP}"
echo -e  "Walltime"

set -o pipefail # important for chained expressions

for sim in $(find $SIMBASE -maxdepth 1 -type d | sort); do
	[[ $sim == "$SIMBASE" ]] && continue
	shortsimname="$(basename "$sim")"
	printf "%-40s" $shortsimname
	if ls $sim/output/*.asc > $h 2>/dev/null; then
		echo -ne "$SEP" $(ls $(<$h) | wc -l) # number of files
		echo -ne "$SEP" $(wc -l $(<$h) | head -n-1 | awk '{ print $1 }' | uniq) # entries in each file
		echo -ne "$SEP" $(wc -l $sim/*log | head -n1 | awk '{ print $1 }') # number of lines in logfile

		if T=$(grep -oh "Finished ExaHyPE successfully" $sim/*.log); then
			echo -ne "$SEP FINISHED"
		elif T=$(grep -oh "ExaHyPE binary failed" $sim/*.log); then
			echo -ne "$SEP FAILED"
		else
			echo -ne "$SEP RUNNING"
		fi

		# filter lines like
		# ... exahype::runners::Runner::startNewTimeStep(...)         step 499 t_min          =0.511616 ...
		if LASTSTEP=$(grep -h startNewTimeStep $sim/*.log | grep step); then
			# red black box, thats why we log for, isnt it?
			BLACKBOX="$(echo "$LASTSTEP" | tail -n1 | awk '{print $6" "$7" "$8" "$9 }')"
		else
			# apparently not even onetime step started.
			BLACKBOX="None started"
		fi
		echo -ne "$SEP $BLACKBOX"

		# `time` from bash normally outputs a line "user	0m0.004s"
		if T=$(grep -E "^user\s+[0-9]" $sim/*.log | awk '{ print $4 }'); then
			echo -ne "$SEP $T"
		# sometimes, we also have something like "321.27user 0.28system 5:23.36elapsed" which is in SECONDS
		elif T=$(grep -Eo "[0-9.]+user\s+" $sim/*.log | tr -d 'user'); then
			T=$(bc <<< "$T / 60" ) # convert seconds to minutes, result is an integer
			echo -ne "$SEP ${T}m" # similarize format
		else
			echo -ne "$SEP noWalltime"
		fi
	else
		echo -n "noResults"
	fi
	echo # newline 
done 
