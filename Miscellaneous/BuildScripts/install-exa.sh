#!/bin/bash

# optional command line argument: how to call the binary
exaname="${1:-exa}"

cd "$(dirname "$0")"


if exaloc=$(which $exaname 2>&1); then
	echo "You already have 'exa' installed at $exaloc"
	ls -l $exaloc
	echo "Use '$0 [name of new exa]' to install exa with an alias (example: 'exa-disk2')"
else
	mkdir -p ~/bin
	buildscripts=$PWD
	cd ~/bin
	ln -s $buildscripts/exa.sh
	mv exa.sh $exaname
	echo "Installed $exaname to ~/bin:"
	ls -l $exaname
	echo "Just use it by typing:"
	which $exaname || { echo "Please log out and log in again to make it work."; }
fi
