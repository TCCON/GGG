
always_yes=false
pyargs=""

for arg in $@; do
    case $arg in 
        -y|--yes)
            always_yes=true
            pyargs="$pyargs --yes"
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
   if [[ $req == 'Y' || $req == 'y' || $always_yes ]] ; then
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
    conda activate ggg-tccon-default
fi

chmod u+x clone_netcdf_writer.sh
./clone_netcdf_writer.sh $pyargs

exit $?
