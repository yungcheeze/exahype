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
		
		# determine logfile
		simlog="$(ls $sim/*.log)"
		if (( "$(ls -1 $sim/*.log | wc -l)" > 1 )); then
			echo -ne "$SEP 'Error, multiple logfiles'"
			>&2 echo "Multiple logfiles for $sim: $simlog"
		else
			if T="$(grep -oh "Finished ExaHyPE successfully" $simlog)"; then
				echo -ne "$SEP FINISHED"
			elif T="$(grep -oh "ExaHyPE binary failed" $simlog)"; then
				echo -ne "$SEP FAILED"
			elif T="$(grep -oh "Peano terminates successfully" $simlog)"; then
				# this can occur when manually run the simulation
				echo -ne "$SEP PEANOFINISHED"
			elif T="$(grep -iE "Signal: Aborted|Asertion .* failed" $simlog)"; then
				# crash due to asserts, memory violations, etc.
				echo -ne "$SEP CRASHED"
			else
				echo -ne "$SEP RUNNING?"
			fi

			# in assertion mode, filter lines like
			# ... exahype::runners::Runner::startNewTimeStep(...)         step 499 t_min          =0.511616 ...
			if LASTSTEP="$(grep -h startNewTimeStep $simlog | grep step)"; then
				# read black box, thats why we log for, isnt it?
				BLACKBOX="$(echo "$LASTSTEP" | tail -n1 | awk '{print $6" "$7" "$8" "$9 }')"
			# in release mode, filter lines like
			# 2017-02-01 12:46:44  158753       info         step 6482        t_min          =3.00505
			elif LASTSTEP="$(grep -E 'info\s+step\s+[0-9]' $simlog | grep t_min)"; then
				BLACKBOX="$(echo "$LASTSTEP" | tail -n1 | awk '{print $5" "$6" "$7" "$8 }')"
			else
				# apparently not even onetime step started.
				BLACKBOX="None started"
			fi
			echo -ne "$SEP $BLACKBOX"

			# `time` from bash normally outputs a line "user	0m0.004s"
			if T="$(grep -E "^user\s+[0-9]" $simlog | awk '{ print $4 }')"; then
				echo -ne "$SEP $T"
			# sometimes, we also have something like "321.27user 0.28system 5:23.36elapsed" which is in SECONDS
			elif T="$(grep -Eo "[0-9.]+user\s+" $simlog | tr -d 'user')"; then
				T="$(bc <<< "$T / 60" )" # convert seconds to minutes, result is an integer
				echo -ne "$SEP ${T}m" # similarize format
			else
				echo -ne "$SEP noWalltime"
			fi
		fi # multiple logfiles
	else
		echo -n "noResults"
	fi # ascii reductions available
	echo # newline 
done 
