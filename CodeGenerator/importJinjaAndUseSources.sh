#!/bin/bash
#
# Use this script to import jinja2 source code locally
# and set the Codegenerator to use it

GIT_URL="https://github.com/pallets/jinja.git"
JINJA_LOCAL_DIR="jinja"

if [ -d "$JINJA_LOCAL_DIR" ] && [ -e "$JINJA_LOCAL_DIR"/.git ]; then
  echo "jinja already imported. Updating..."
  cd "$JINJA_LOCAL_DIR"
  git pull
  cd ..
else 
  echo "Cloning jinja2 from $GIT_URL"
  git clone "$GIT_URL" "$JINJA_LOCAL_DIR"
fi

echo "Set TemplatingUtils.py to use the local jinja sources"
sed -i -e "s#isJinjaAvailableAsPackage=True#isJinjaAvailableAsPackage=False#" TemplatingUtils.py

echo "You may want to stop tracking local changes to TemplatingUtils.py with: 'git update-index --skip-worktree TemplatingUtils.py'"