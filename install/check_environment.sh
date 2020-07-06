#!/bin/bash

# Necessary to use conda commands inside this script, even though this was already 
# done in master.sh.
source $GGGPATH/install/.init_conda

# make this an array of environments so that we can check each one
# to see if it is ggg-tccon-default
conda_envs=($(conda env list | awk '{print $1}'))
env_exists=false
for e in ${conda_envs[*]}; do
    if [[ $e == ggg-tccon-default ]]; then
        env_exists=true
        break
    fi  
done

if $env_exists; then
        echo "The ggg-tccon-default python environment already exists; update from environment.yml"
        conda env update -f environment.yml
else
        echo "The ggg-tccon-default python environment does not exist"
        echo "Creating ggg-tccon-default python environment from environment.yml"
        conda env create -f environment.yml
fi
