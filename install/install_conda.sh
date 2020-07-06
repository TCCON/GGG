#!/bin/bash

if [ $# -gt 0 ]; then
    conda_dest=$1
else
    conda_dest=$HOME
fi

if [ "$(uname)" == "Linux" ] ; then
    conda="Anaconda3-2019.10-Linux-x86_64.sh"
elif [ "$(uname)" == "Darwin" ] ; then
    conda="Anaconda3-2019.10-Linux-x86_64.sh"
else
    echo "$(uname) platform not supported. Please download and install Anaconda3 manually."
    exit 1
fi

cd $conda_dest
echo "Will install $conda in `pwd`."
read -p "Continue? " answer
case $answer in
    [yY])
        echo "Beginning Anaconda3 download and install..."
        ;;
    *)
        echo "Aborting Anaconda3 install."
        echo "To install to a different location, pass that path as the command line argument."
        exit 1
        ;;
    esac


if [ ! -f $conda ]; then
    echo "Downloading $conda"
    wget https://repo.continuum.io/archive/$conda
else
    echo "Found existing $conda, reusing"
fi

echo "Installing $conda"
exit 0
bash $conda
