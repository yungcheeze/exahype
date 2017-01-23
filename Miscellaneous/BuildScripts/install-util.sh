#!/bin/bash

# A small "installation" utility. Usage:
#  install.sh <NameOfSourceFile> <NameOftargetFile>

# mandatory command line argument: name of source file
sourcefile="$1"

# mandatory command line argument: how to call the binary
targetfile="$2"

# where to install
bin="$HOME/bin"

cd "$(dirname "$0")"


if targetloc=$(which $targetfile 2>&1); then
	echo "You already have $sourcefile installed at $targetloc"
	ls -l $targetloc
	#echo "Use '$0 [name of new exa]' to install exa with an alias (example: 'exa-disk2')"
	exit -1
else
	mkdir -p $bin && cd $bin
	ln -s $sourcefile $targetfile
	echo "Installed $sourcefile to $bin:"
	ls -l $targetfile
	echo "Just use it by typing:"
	which $targetfile || { echo "Please log out and log in again to make it work."; }
fi
