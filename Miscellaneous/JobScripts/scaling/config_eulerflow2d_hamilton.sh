#!/bin/bash
export SCHEDULER_INTERACTIVE=1
# exahype shared memory properties
export SCHEDULER_PROJECT_DIR
# relative from working dir or absolute
export SCHEDULER_PROJECT_DIR=eulerflow2d
export SCHEDULER_SCRIPT=./scheduler_submit_job_hamilton.sh
export SCHEDULER_OUTPUT_DIR=${SCHEDULER_PROJECT_DIR}/scaling
# relative from project dir
export SCHEDULER_EXECUTABLE=./ExaHyPE-Euler2d-TBB
export SCHEDULER_EXECUTABLE_SERIAL=./ExaHyPE-Euler2d-None
export SCHEDULER_SPEC_FILE=../eulerflow2d.exahype
