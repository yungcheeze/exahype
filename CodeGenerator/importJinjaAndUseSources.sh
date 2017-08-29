#!/bin/bash
#
# Use this script to import jinja2 source code locally
# and set the Codegenerator to use it

# Configuration variable
GIT_URL="https://github.com/pallets/jinja.git"
JINJA_LOCAL_DIR="jinja"

# local var to resolve relative path correctly
scriptDir=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
currentLocation=$(pwd)

# move to the CodeGenerator directory
cd $scriptDir

# import or update jinja2 to CodeGenerator/jinja
if [ -d "$JINJA_LOCAL_DIR" ] && [ -e "$JINJA_LOCAL_DIR"/.git ]; then
  echo "jinja already imported. Updating..."
  cd "$JINJA_LOCAL_DIR"
  #git pull
  cd ..
else 
  echo "Cloning jinja2 from $GIT_URL"
  #git clone "$GIT_URL" "$JINJA_LOCAL_DIR"
fi

# modify the Codegenerator default configuration
echo "Set TemplatingUtils.py to use the local jinja sources"
sed -i -e "s#isJinjaAvailableAsPackage=True#isJinjaAvailableAsPackage=False#" TemplatingUtils.py

echo "You may want to stop tracking local changes to TemplatingUtils.py with: 'git update-index --skip-worktree TemplatingUtils.py'"

# move back to where the script was called
cd $currentLocation
