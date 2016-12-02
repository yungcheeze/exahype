#!/bin/bash
#
# Use this script as a loose coupling of the ExaHyPE-Engine git
# repository to the Peano svn repository.
#

SVN_URL="svn://svn.code.sf.net/p/peano/code/trunk/src"

if test -e .svn; then
	echo "Updating peano"
	svn update
else
	echo "Checking out peano from $SVN_URL"
	svn checkout $SVN_URL .
fi
