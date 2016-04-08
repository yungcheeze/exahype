#!/bin/bash
# General config file. This file is included by the specialised config files.
export SCHEDULER_INTERACTIVE=1
# relative from working dir or absolute
export SCHEDULER_PROJECT_DIR=eulerflow2d
export SCHEDULER_SCRIPT=./scheduler_submit_job.sh
# relative from project dir
export SCHEDULER_EXECUTABLE=./ExaHyPE-Euler2d-TBB
export SCHEDULER_EXECUTABLE_SERIAL=./ExaHyPE-Euler2d-None
export SCHEDULER_OUTPUT_DIR=${SCHEDULER_PROJECT_DIR}/scaling

export SCHEDULER_OUTPUT_PREFIX=1tr_dummy
export SCHEDULER_SPEC_FILE=../eulerflow2d.exahype
export SHAREDMEM=TBB
BASH_ENV=params_sandybridge.sh ./scheduler.sh
export SHAREDMEM=None
BASH_ENV=params_sandybridge.sh ./scheduler_serial.sh

export SCHEDULER_OUTPUT_PREFIX=4tr_dummy
export SCHEDULER_SPEC_FILE=../eulerflow2d_4tr.exahype
export SHAREDMEM=TBB
BASH_ENV=params_sandybridge.sh ./scheduler.sh
export SHAREDMEM=None
BASH_ENV=params_sandybridge.sh ./scheduler_serial.sh

