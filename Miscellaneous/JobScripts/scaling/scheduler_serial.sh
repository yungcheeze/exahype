#!/bin/bash
###############################################################################
# JOB SCHEDULER
# Writes core counts in ExaHyPE project specification file (*.exahype) and
# runs an ExaHyPE project executable.
###############################################################################
# User needs to provide the following parameters
# CC                          (string)
# MODE                        (string)
# SHAREDMEM                   (string)
# SCHEDULER_PROJECT_DIR       (string)
# SCHEDULER_SCRIPT            (string)
# SCHEDULER_EXECUTABLE_SERIAL (string)
# SCHEDULER_SPEC_FILE         (string)
# SCHEDULER_PROCESSOR_ID      (string)
# SCHEDULER_OUTPUT_PREFIX     (string)
# SCHEDULER_OUTPUT_DIR        (string)
# SCHEDULER_MPI               (array)
# SCHEDULER_CORES             (array)
# SCHEDULER_RUNS              (number)

###############################################################################
# READ IN PARAMETERS
###############################################################################
printf "###############################################################################\n"
printf "SCRIPT STARTED WITH THE FOLLOWING PARAMETERS:\n"
printf "###############################################################################\n"
printf "CC=${CC}\n"
printf "MODE=${MODE}\n"
printf "SHAREDMEM=${SHAREDMEM}\n\n"
printf "SCHEDULER_PROJECT_DIR=${SCHEDULER_PROJECT_DIR}\n"
printf "SCHEDULER_SCRIPT=${SCHEDULER_SCRIPT}\n"
printf "SCHEDULER_EXECUTABLE_SERIAL=${SCHEDULER_EXECUTABLE_SERIAL}\n"
printf "SCHEDULER_SPEC_FILE=${SCHEDULER_SPEC_FILE}\n"
printf "SCHEDULER_SYSTEM_ID=${SCHEDULER_SYSTEM_ID}\n"
printf "SCHEDULER_OUTPUT_PREFIX=${SCHEDULER_OUTPUT_PREFIX}\n"
printf "SCHEDULER_OUTPUT_DIR=${SCHEDULER_OUTPUT_DIR}\n"
printf "SCHEDULER_RUNS=${SCHEDULER_RUNS}\n"
printf "SCHEDULER_INTERACTIVE=${SCHEDULER_INTERACTIVE}\n\n"
printf "\n"

if [ "${SHAREDMEM}" != "None" ]
then
  printf "ERROR: SHAREDMEM must be set to \'None\'!\n"
  printf "       Please make also sure your executable is compiled with\n"
  printf "       this option!\n"
  printf "       Script execution failed!\n"
  exit 1
fi

printf "Turn interactive mode on or off by setting SCHEDULER_INTERACTIVE=\'1\' or \n\'0\', respectively.\n"
if [ "${SCHEDULER_INTERACTIVE}" == "1" ]
then
printf "###############################################################################\n"
  printf "NOTE: Check the parameters make sure that they are correct\n      and not empty!\n\n"
  read -p "Press [ Ctrl+C ] to exit or any other key to continue!"
fi

###############################################################################
# CREATE OUTPUT DIRECTORIES
###############################################################################
PREFIX=$(date +"%y%m%d")_${SCHEDULER_OUTPUT_PREFIX}
OUTPUT_DIR=${SCHEDULER_OUTPUT_DIR}/${PREFIX}

printf "Trying to create output dir \'${SCHEDULER_OUTPUT_DIR}\' \n"
mkdir ${SCHEDULER_OUTPUT_DIR}

printf "Trying to create prefixed subfolder \'${OUTPUT_DIR}\' \n"
mkdir ${OUTPUT_DIR}

###############################################################################
# RUN THE EXPERIMENTS
###############################################################################
for t in 1
do
  # Create temporary spec file with new core count

  ( \
  cd ${SCHEDULER_PROJECT_DIR} && \
  sed "s/cores                    =.*/cores                    = ${t}/" ${SCHEDULER_SPEC_FILE} > ${SCHEDULER_SPEC_FILE}_$t \
  )
  
for r in $(seq 1 ${SCHEDULER_RUNS})
  do  
    OUTPUT_FILE=${OUTPUT_DIR}/${PREFIX}_n1_t${t}_r${r}_${CC}_${SHAREDMEM}.txt
      
    printf "Writing file \'$OUTPUT_FILE\'... "
    ${SCHEDULER_SCRIPT} ${SCHEDULER_PROJECT_DIR} ${SCHEDULER_EXECUTABLE_SERIAL} ${SCHEDULER_SPEC_FILE}_$t ${OUTPUT_FILE}
    printf "done!\n"
  done
  # Delete temporary spec file 
  #(cd ${SCHEDULER_PROJECT_DIR}  &&  rm ${SCHEDULER_SPEC_FILE}_$t)
done
printf "Script finished successfully!\n"
