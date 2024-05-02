usage () {
    echo "$0 [ -y | --yes ]"
    echo " The -y or --yes flag will automatically answer 'y' to any interactive question."
    echo " This script also respects the environmental variable GGG_USE_MICROMAMBA."
    echo " If GGG_USE_MICROMAMBA=1, then the netCDF installer will be configured to use"
    echo " micromamba instead of conda *if* it is found on your path. To force the use"
    echo " of micromamba even if it does not appear to be on your path, set GGG_USE_MICROMAMBA=2"
    echo " instead."
}
always_yes=false
pyargs=""

for arg in $@; do
    case $arg in 
        -y|--yes)
            always_yes=true
            pyargs="$pyargs --yes"
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
    echo "Could not configure GGG to use Anaconda3, aborting."
    exit 1
fi

chmod u+x check_environment.sh
./check_environment.sh $pyargs

initfile=$GGGPATH/install/.init_conda
if [ ! -f $initfile ]; then
    echo "Could not initialize conda, $initfile does not exist. Has $GGGPATH/install/check_python.sh been run?"
    exit 1
else
    source $initfile
    if [ -z $NO_CONDA_ACTIVATE ] || [ $NO_CONDA_ACTIVATE == 0 ]; then
        echo "Activating $GGGPATH/install/.condaenv"
        $CONDACMD activate "$GGGPATH/install/.condaenv"
        echo $CONDA_PREFIX
    fi
fi

chmod u+x clone_netcdf_writer.sh
./clone_netcdf_writer.sh $pyargs

exit $?
