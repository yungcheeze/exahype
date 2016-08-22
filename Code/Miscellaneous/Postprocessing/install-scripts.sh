#!/bin/sh
#
# A small script to install softlinks to your Home's bin or
# to /usr/local/bin or where you would like to have these
# scripts at.

install_dir=$HOME/bin
Postprocessing=$PWD

ln -vs $Postprocessing/exaplayer.py $install_dir/exaplayer
ln -vs $Postprocessing/exareader.py $install_dir/exareader

# aux tools
ln -vs $Postprocessing/webmize $install_dir/webmize

