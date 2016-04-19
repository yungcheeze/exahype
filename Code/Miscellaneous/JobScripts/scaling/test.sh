#!/bin/bash
source ${EXAHYPE_PARAMS_FILE}

echo $EXAHYPE_SYSTEM_ID

printf "EXAHYPE_MPI_PROCESSES="
for n in "${EXAHYPE_MPI_PROCESSES[@]}"; do printf "$n,"; done
printf "\n"
printf "EXAHYPE_CORES="
for t in "${EXAHYPE_CORES[@]}"; do printf "$t," ; done
printf "\n"

