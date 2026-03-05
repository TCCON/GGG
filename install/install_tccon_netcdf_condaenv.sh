#!/bin/bash

# This script is intended to be run in the tccon_netcdf directory

# In some cases, an active environment does not persist when a subshell
# is created by executing a script so we double check that here.
# Use realpath on both to avoid any comparisons failing because one uses
# the real physical path and the other goes through symlinks.

if [ -z $NO_CONDA_ACTIVATE ] || [ $NO_CONDA_ACTIVATE == 0 ]; then
  conda_prefix=$(realpath $CONDA_PREFIX)
  expected_prefix=$(realpath "${GGGPATH}/install/.condaenv")
  if [ $conda_prefix != $expected_prefix ] ; then
    echo "The $expected_prefix python environment is not active. CONDA_PREFIX=$conda_prefix. Cannot proceed." 
    exit 1
  fi
fi

# --editable: do not copy to the environment, just run from here
# --isolated: ignore any user configuration, which can sometimes
#             redirect the installation to weird places.
# --no-deps: do not install dependencies with pip. They should all
#            be managed by the conda environment.
$PIPCMD install --no-deps --isolated --editable .
$PYTHONCMD copy_scripts.py --env-prefix "${GGGPATH}/install/.condaenv/" "${GGGPATH}/bin/"
