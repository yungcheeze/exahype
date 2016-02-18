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
printf "EXAHYPE_RUNS=${EXAHYPE_RUNS}\n"
printf "\n"

if [ "${SHAREDMEM}" != "None" ]
then
  printf "ERROR: SHAREDMEM must be set to \'None\'!\n"
  printf "       Please make also sure your executable is compiled with\n"
  printf "       this option!\n"
  printf "       Script execution failed!\n"
  exit 1
fi

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
for t in 1
do
  # Create temporary spec file with new core count
  sed "s/cores                    =.*/cores                    = ${t}/" ${EXAHYPE_SPEC_FILE} > ${EXAHYPE_SPEC_FILE}_$t
  for r in $(seq 1 ${EXAHYPE_RUNS})
  do  
    OUTPUT_FILE=${OUTPUT_DIR}/${PREFIX}_n1_t${t}_r${r}_${CC}_${SHAREDMEM}.txt
      
    printf "Writing file \'$OUTPUT_FILE\'... "
    # RUN THE EXECUTABLE AND PIPE IT INTO OUTPUT FILE 
    # @VASCO: YOU HAVE TO APPLY YOUR CHANGES HERE!
    ${EXAHYPE_EXECUTABLE} ${EXAHYPE_SPEC_FILE}_$t > ${OUTPUT_FILE}
    printf "done!\n"
  done
  # Delete temporary spec file 
  rm ${EXAHYPE_SPEC_FILE}_$t
done
printf "Script finished successfully!\n"
