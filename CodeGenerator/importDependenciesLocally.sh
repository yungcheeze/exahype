#!/bin/bash
#
# Use this script to import jinja2 source code locally (with its dependency)

# Configuration variable
JINJA_GIT_URL="https://github.com/pallets/jinja.git"
MARKUPSAFE_GIT_URL="https://github.com/pallets/markupsafe.git"
JINJA_LOCAL_DIR="dependencies/jinja"
MARKUPSAFE_LOCAL_DIR="dependencies/markupsafe-install"

# local var to resolve relative path correctly
scriptDir=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
currentLocation=$(pwd)

# move to the CodeGenerator directory
cd "$scriptDir"

# import or update jinja2 to CodeGenerator/dependencies/jinja
if [ -d "$JINJA_LOCAL_DIR" ] && [ -e "$JINJA_LOCAL_DIR"/.git ]; then
  echo "jinja already imported. Updating..."
  cd "$JINJA_LOCAL_DIR"
  git pull
  cd "$scriptDir"
else 
  echo "Cloning jinja2 from $JINJA_GIT_URL"
  git clone "$JINJA_GIT_URL" "$JINJA_LOCAL_DIR"
fi

# import or update markup-safe to CodeGenerator/dependencies/markup-safe
if [ -d "$MARKUPSAFE_LOCAL_DIR" ] && [ -e "$MARKUPSAFE_LOCAL_DIR"/.git ]; then
  echo "jinja already imported. Updating..."
  cd "$MARKUPSAFE_LOCAL_DIR"
  git pull
  cd "$scriptDir"
else 
  echo "Cloning markup-safe from $MARKUPSAFE_GIT_URL"
  git clone "$MARKUPSAFE_GIT_URL" "$MARKUPSAFE_LOCAL_DIR"
fi

# make sym links
ln -s "$MARKUPSAFE_LOCAL_DIR"/markupsafe
ln -s "$JINJA_LOCAL_DIR"/jinja2

# move back to where the script was called
cd "$currentLocation"
