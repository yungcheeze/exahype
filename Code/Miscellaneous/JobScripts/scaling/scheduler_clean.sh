#!/bin/bash
###############################################################################
# JOB SCHEDULER
# Delete temporary spec files.
###############################################################################
# User needs to provide the following parameters
# SCHEDULER_PROJECT_DIR   (string)
# SCHEDULER_SPEC_FILE     (string)
# SCHEDULER_CORES         (array)
for t in "${SCHEDULER_CORES[@]}"
do
  (cd ${SCHEDULER_PROJECT_DIR}  &&  rm ${SCHEDULER_SPEC_FILE}_$t)    
done
