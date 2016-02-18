#!/bin/bash
###############################################################################
# JOB SCHEDULER
# Writes core counts in ExaHyPE project specification file (*.exahype) and
# runs an ExaHyPE project executable.
###############################################################################
# User needs to provide the following parameters
# CC                    (string)
# MODE                  (string)
# SHAREDMEM             (string)
# EXAHYPE_PROJECT_DIR   (string)
# EXAHYPE_EXECUTABLE    (string)
# EXAHYPE_SPEC_FILE     (string)
# EXAHYPE_PROCESSOR_ID  (string)
# EXAHYPE_OUTPUT_PREFIX (string)
# EXAHYPE_OUTPUT_DIR    (string)
# EXAHYPE_MPI           (array)
# EXAHYPE_CORES         (array)
# EXAHYPE_NRUNS         (number)

###############################################################################
# READ IN PARAMETERS
###############################################################################
printf "###############################################################################\n"
printf "SCRIPT STARTED WITH THE FOLLOWING PARAMETERS:\n"
printf "###############################################################################\n"
printf "CC=${CC}\n"
printf "MODE=${MODE}\n"
printf "SHAREDMEM=${SHAREDMEM}\n\n"
printf "EXAHYPE_PROJECT_DIR=${EXAHYPE_PROJECT_DIR}\n"
printf "EXAHYPE_EXECUTABLE=${EXAHYPE_EXECUTABLE}\n"
printf "EXAHYPE_SPEC_FILE=${EXAHYPE_SPEC_FILE}\n"
printf "EXAHYPE_SYSTEM_ID=${EXAHYPE_SYSTEM_ID}\n"
printf "EXAHYPE_PREFIX=${EXAHYPE_PREFIX}\n"
printf "EXAHYPE_OUTPUT_DIR=${EXAHYPE_OUTPUT_DIR}\n"
printf "EXAHYPE_INTERACTIVE=${EXAHYPE_INTERACTIVE}\n\n"

# Select MPI process and core counts based on system id.
# NOTE: I can't outsource these numbers since bash does not support exporting
# arrays out of the box.
if [ "${EXAHYPE_SYSTEM_ID}" == "sandybridge" ]
then
  EXAHYPE_MPI_PROCESSES=(1)
  EXAHYPE_CORES=(2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 32)
elif [ "${EXAHYPE_SYSTEM_ID}" == "haswell" ]
then
  EXAHYPE_MPI_PROCESSES=(1)
  EXAHYPE_CORES=(2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 56)
elif [ "${EXAHYPE_SYSTEM_ID}" == "mic" ]
then
  export EXAHYPE_MPI_PROCESSES=(1)
  export EXAHYPE_CORES=(1 2 8 12 16 18 24 60 72 120 180 240)
else
  printf "System \'${EXAHYPE_SYSTEM_ID}\' is not supported!"
  exit 1;
fi

printf "EXAHYPE_MPI_PROCESSES="
for n in "${EXAHYPE_MPI_PROCESSES[@]}"; do printf "$n,"; done
printf "\n"
printf "EXAHYPE_CORES="
for t in "${EXAHYPE_CORES[@]}"; do printf "$t," ; done
printf "\n"
printf "EXAHYPE_RUNS=${EXAHYPE_RUNS}\n"
printf "\n"

printf "Turn interactive mode on or off by setting EXAHYPE_INTERACTIVE=\'1\' or \n\'0\', respectively.\n"
if [ "${EXAHYPE_INTERACTIVE}" == "1" ]
then
printf "###############################################################################\n"
  printf "NOTE: Check the parameters make sure that they are correct\n      and not empty!\n\n"
  read -p "Press [ Ctrl+C ] to exit or any other key to continue!"
fi

###############################################################################
# CREATE OUTPUT DIRECTORIES
###############################################################################
PREFIX=$(date +"%y%m%d")_${EXAHYPE_PREFIX}
OUTPUT_DIR=${EXAHYPE_OUTPUT_DIR}/${PREFIX}

printf "Trying to create scaling subfolder \'${EXAHYPE_OUTPUT_DIR}\' \n"
mkdir ${EXAHYPE_OUTPUT_DIR}

printf "Trying to create output directory \'${OUTPUT_DIR}\' \n"
mkdir ${OUTPUT_DIR}

###############################################################################
# RUN THE EXPERIMENTS
###############################################################################
for n in "${EXAHYPE_MPI_PROCESSES[@]}";
do
  for t in "${EXAHYPE_CORES[@]}"
  do
    # Create temporary spec file with new core count
    sed "s/cores                    =.*/cores                    = ${t}/" ${EXAHYPE_SPEC_FILE} > ${EXAHYPE_SPEC_FILE}_$t
    for r in $(seq 1 ${EXAHYPE_RUNS})
    do  
      OUTPUT_FILE=${OUTPUT_DIR}/${PREFIX}_n${n}_t${t}_r${r}_${CC}_${SHAREDMEM}.txt
      
      printf "Writing file \'$OUTPUT_FILE\'... "
      # RUN THE EXECUTABLE AND PIPE IT INTO OUTPUT FILE 
      # @VASCO: YOU HAVE TO APPLY YOUR CHANGES HERE!
      ${EXAHYPE_EXECUTABLE} ${EXAHYPE_SPEC_FILE}_$t > ${OUTPUT_FILE}
      printf "done!\n"
    done
    # Delete temporary spec file 
    rm ${EXAHYPE_SPEC_FILE}_$t
  done
done
printf "Script finished successfully!\n"
