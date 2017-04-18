#!/bin/bash

# optional command line argument: how to call the binary
exaname="${1:-exa}"

# in order to install at another location, use like
#   target="/usr/local/bin" .../install-exa.sh"
target="${target:=$HOME/bin}"

cd "$(dirname "$0")"


if exaloc=$(which $exaname 2>&1); then
	echo "You already have 'exa' installed at $exaloc"
	ls -l $exaloc
	echo "Use '$0 [name of new exa]' to install exa with an alias (example: 'exa-disk2')"
else
	mkdir -p $target
	buildscripts=$PWD
	cd $target
	ln -s $buildscripts/exa.sh
	mv exa.sh $exaname
	echo "Installed $exaname to $target:"
	ls -l $exaname
	echo "Just use it by typing:"
	which $exaname || { echo "Please log out and log in again to make it work."; }
fi
