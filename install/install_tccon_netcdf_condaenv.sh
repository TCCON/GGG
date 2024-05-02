#!/bin/bash

# In some cases, an active environment does not persist when a subshell
# is created by executing a script so we double check that here.
# Use realpath on both to avoid any comparisons failing because one uses
# the real physical path and the other goes through symlinks.

if [ -z $NO_CONDA_ACTIVATE ] || [ $NO_CONDA_ACTIVATE == 0 ]; then
  conda_prefix=$(realpath $CONDA_PREFIX)
  expected_prefix=$(realpath $GGGPATH/install/.condaenv)
  if [ $conda_prefix != $expected_prefix ] ; then
    echo "The $expected_prefix python environment is not active. CONDA_PREFIX=$conda_prefix. Cannot proceed." 
    exit 1
  fi
fi

if [ ! -d scripts ]; then
    mkdir -v scripts
fi

# develop: do not copy to $SITEDIR/, just run from here
# --no-user-cfg: ignore any ~/.pydistutils.cfg file
# --script-dir: write command line scripts to the given directory
$PYTHONCMD setup.py --no-user-cfg develop --script-dir=./scripts

if [ $? == 0 ]; then
    [ -f ./scripts/write_netcdf ] && mv -v ./scripts/write_netcdf $GGGPATH/bin/
    [ -f ./scripts/compare_netcdf ] && mv -v ./scripts/compare_netcdf $GGGPATH/bin/
    [ -f ./scripts/concat_netcdf ] && mv -v ./scripts/concat_netcdf $GGGPATH/bin/
    [ -f ./scripts/subset_netcdf ] && mv -v ./scripts/subset_netcdf $GGGPATH/bin/
    [ -f ./scripts/update_site_info ] && mv -v ./scripts/update_site_info $GGGPATH/bin/
    [ -f ./scripts/update_manual_flags ] && mv -v ./scripts/update_manual_flags $GGGPATH/bin/
    exit 0
else
    echo "setup.py for tccon_netcdf failed"
    exit 1
fi
