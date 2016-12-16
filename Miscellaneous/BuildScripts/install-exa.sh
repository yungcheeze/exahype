#!/bin/bash

if exaloc=$(which exa 2>&1); then
	echo "You already have 'exa' installed at $exaloc"
	ls -l $exaloc
	echo "Don't do anything."
	# todo: Could install another exa with an appendix
else
	mkdir -p ~/bin
	buildscripts=$PWD
	cd ~/bin
	ln -s $buildscripts/exa.sh
	mv exa.sh exa
	echo "Installed exa to ~/bin:"
	ls -l exa
	echo "Just use it by typing:"
	which exa || { echo "Please log out and log in again to make it work."; }
fi
