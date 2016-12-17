#!/bin/bash

verbose () { echo $@; $@; }

[[ -e simulations ]] && verbose cd "simulations"

# alternatively, use find:
#   find . -iname *.env -or -iname *.asc -or -iname *.log
# it's slower on my system than globbing

verbose tar cfz simulation-tables.tar.gz */output/ */*.{env,log} \
	&& echo "done" || echo "failure"
