#!/bin/bash -x

# this deletes the crazy amount of log files created by ExaHyPE
# which makes it impossible even to list directories.

find -L . -iname *.log-file -exec rm {} \;
