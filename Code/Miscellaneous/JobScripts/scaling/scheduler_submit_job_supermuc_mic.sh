#!/bin/bash
###############################################################################
# Submits a job.
###############################################################################
# Is invoked with the following arguments
# $1 absolute path to the folder containing the executable
# $2 the executable (./)
# $3 relative path to the specification file from viewpoint of executable
#    or absolute path
# $4 path to the output file

JOB_SCRIPT=scheduler_job_supermuc_mic.job

sed "s,{PROJECT_DIR},$1,g" ${JOB_SCRIPT} > ${JOB_SCRIPT}_tmp
sed "s,{RUN},$2 $3,g" ${JOB_SCRIPT}_tmp > ${JOB_SCRIPT}_tmp1
sed "s,{OUTPUT_FILE},$4,g" ${JOB_SCRIPT}_tmp1 > ${JOB_SCRIPT}_tmp

rm ${JOB_SCRIPT}_tmp1

sbatch ${JOB_SCRIPT}_tmp

rm ${JOB_SCRIPT}_tmp
