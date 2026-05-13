#!/bin/bash
usage () {
    echo "$0 [ -y | --yes ] [--allow-env-mismatch]"
    echo " The -y or --yes flag will automatically answer 'y' to any interactive question."
    echo ""
    echo " The --allow-env-mismatch flag allows the installation to continue even if the "
    echo " checksum of the environment.yml file does not match the last time this script "
    echo " was run. This will not identify manual changes to the environment, only changes "
    echo " across GGG versions."
    echo ""
    echo " This script also respects the environmental variables GGG_USE_MICROMAMBA and GGG_USE_PIP."
    echo " If GGG_USE_MICROMAMBA=1, then the netCDF installer will be configured to use"
    echo " micromamba instead of conda *if* it is found on your path. To force the use"
    echo " of micromamba even if it does not appear to be on your path, set GGG_USE_MICROMAMBA=2"
    echo " instead. If GGG_USE_PIP=1, then this will use standard Python virtual environments."
    echo " See https://tccon-wiki.caltech.edu/Main/InstallingPython for more information."
}
always_yes=false
enforce_env_match=true
pyargs=""

for arg in $@; do
    case $arg in 
        -y|--yes)
            always_yes=true
            pyargs="$pyargs --yes"
            ;;
        --allow-env-mismatch)
            enforce_env_match=false
            ;;
        -h|--help)
            usage
            exit 0
            ;;
    esac
done

echo " Installing GGG"
echo " Using GGGPATH =" $GGGPATH
echo "$GGGPATH/install == $(pwd)"
if [ ! $GGGPATH/install == `pwd` ] ; then
   echo " Your current directory: `pwd`"
   echo " Does not match the GGGPATH install directory."
   read -p " Continue? (Y/N) " req
   if [[ $req == 'Y' ]] || [[ $req == 'y' ]] || $always_yes ; then
      echo " Continuing..."
   else
      echo " Quitting install. Please change your GGGPATH."
      exit
   fi
fi

chmod u+x check_python.sh
./check_python.sh $pyargs
pyexit=$?
if [ $pyexit != 0 ] ; then
    echo "Could not configure Python for GGG to use, aborting."
    exit 1
fi

chmod u+x check_environment.sh
./check_environment.sh $pyargs
envexit=$?
if [ $envexit != 0 ]; then
    if $enforce_env_match; then
        echo "ERROR: Python environment out of date, aborting netCDF writer install"
        exit $envexit
    else
        echo "WARNING: Python environment likely out of date, but proceeding with Python installation anyway"
    fi
fi

initfile=$GGGPATH/install/.init_conda
if [ ! -f $initfile ]; then
    echo "Could not initialize conda, $initfile does not exist. Has $GGGPATH/install/check_python.sh been run?"
    exit 1
else
    source $initfile
    if [ -z $NO_CONDA_ACTIVATE ] || [ $NO_CONDA_ACTIVATE == 0 ]; then
        $CONDACMD activate "$GGGPATH/install/.condaenv"
    fi
fi

chmod u+x clone_netcdf_writer.sh
./clone_netcdf_writer.sh $pyargs

exit $?
