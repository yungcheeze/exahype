#!/bin/bash
#
# Call this similarly as showSimulationProgress.sh,
# ie. with the SIMBASE like in

SIMBASE="${SIMBASE:=simulations/}"
LIBCONV="$(dirname "$0")"
EXTRACT_DTMIN="$LIBCONV/extract_dtmin.py"

# todo: This should be read consistently out of an simulation parameters.env file.
TIMESTEPFILENAME="timesteps.csv"

# ie. out from the directory where the simulations live in.

for sim in $(find $SIMBASE -maxdepth 1 -type d | sort); do
	[[ $sim == "$SIMBASE" ]] && continue
	shortsimname="$(basename "$sim")"
	printf "%-40s" $shortsimname

	# determine logfile
	simlog="$(ls $sim/*.log)"
	if (( "$(ls -1 $sim/*.log | wc -l)" > 1 )); then
		echo -ne "Error, multiple logfiles"
		>&2 echo "Multiple logfiles for $sim: $simlog"
	else
		echo "Reading out timesteps..."
		$EXTRACT_DTMIN $simlog > $sim/$TIMESTEPFILENAME &
	fi
done

wait
