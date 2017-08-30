#!/bin/bash
#
# Use this script to import jinja2 source code locally (with its dependency)



#************************************
#***** Configuration variables ******
#************************************

#Change this if Jinja2 is available with your python3.
JINJA2_ALREADY_AVAILABLE=false #false => install Jinja2 and MarkupSafe locally, true => skip it. Default value is false.
#git paths
JINJA_GIT_URL="https://github.com/pallets/jinja.git"
MARKUPSAFE_GIT_URL="https://github.com/pallets/markupsafe.git"
LIBXSMM_GIT_URL="https://github.com/hfp/libxsmm.git"
#local clone paths
JINJA_LOCAL_DIR="dependencies/jinja"
MARKUPSAFE_LOCAL_DIR="dependencies/markupsafe"
LIBXSMM_LOCAL_DIR="dependencies/libxsmm"




#************************************
#******* Installation steps *********
#************************************

# local var to resolve relative path correctly
scriptDir=$(dirname -- "$(readlink -f -- "$BASH_SOURCE")")
currentLocation=$(pwd)

# move to the CodeGenerator directory
cd "$scriptDir"

if [ "$JINJA2_ALREADY_AVAILABLE" = false ] ; then
  echo "Importing Jinja2 and MarkupSafe"
  # import or update jinja2 to CodeGenerator/dependencies/jinja
  if [ -d "$JINJA_LOCAL_DIR" ] && [ -e "$JINJA_LOCAL_DIR"/.git ]; then
    echo "Jinja2 already imported. Updating..."
    cd "$JINJA_LOCAL_DIR"
    git pull
    cd "$scriptDir"
  else 
    echo "Cloning Jinja2 from $JINJA_GIT_URL"
    git clone "$JINJA_GIT_URL" "$JINJA_LOCAL_DIR"
  fi

  # import or update markupsafe to CodeGenerator/dependencies/markupsafe
  if [ -d "$MARKUPSAFE_LOCAL_DIR" ] && [ -e "$MARKUPSAFE_LOCAL_DIR"/.git ]; then
    echo "MarkupSafe already imported. Updating..."
    cd "$MARKUPSAFE_LOCAL_DIR"
    git pull
    cd "$scriptDir"
  else 
    echo "Cloning MarkupSafe from $MARKUPSAFE_GIT_URL"
    git clone "$MARKUPSAFE_GIT_URL" "$MARKUPSAFE_LOCAL_DIR"
  fi
else
  echo "Jinja2 is already available with python3"
fi

# import or update libxsmm to CodeGenerator/dependencies/libxsmm
if [ -d "$LIBXSMM_GIT_URL" ] && [ -e "$LIBXSMM_GIT_URL"/.git ]; then
  echo "LIBXSMM already imported. Updating..."
  cd "$LIBXSMM_LOCAL_DIR"
  git pull
  cd "$scriptDir"
else 
  echo "Cloning LIBXSMM from $LIBXSMM_GIT_URL"
  git clone "$LIBXSMM_GIT_URL" "$LIBXSMM_LOCAL_DIR"
fi

# build libxsmm
echo "Build libxsmm gemm generator"
cd "$LIBXSMM_LOCAL_DIR"
make realclean
make generator
cd "$scriptDir"

# make sym links
echo "Make symbolic links"
if [ "$JINJA2_ALREADY_AVAILABLE" = false ] ; then
  ln -s "$MARKUPSAFE_LOCAL_DIR"/markupsafe
  ln -s "$JINJA_LOCAL_DIR"/jinja2
fi
ln -s "$LIBXSMM_LOCAL_DIR"/bin/libxsmm_gemm_generator

# move back to where the script was called
cd "$currentLocation"
