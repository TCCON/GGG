#!/bin/bash

usage() {
    echo "$0 [ -y | --yes ] [ --no-py ] [ NPROCS ] "
    echo "  -y | --yes : Assume user enters YES for all prompts."
    echo "               WARNING: this effectively disables the"
    echo "               GGGPATH check."
    echo "  --no-py : Do not reinstall the Python packages needed, "
    echo "            only recompile/retest GGG."
}


##### Command line argument parsing ####
# Defaults first
nprocessors=1
install_py=true
always_yes=false
pyargs=""

# Now parse
for arg in $@; do
    case $arg in 
        -h|--help)
            usage
            exit 0
            ;;
        -y|--yes)
            always_yes=true
            pyargs="$pyargs --yes"
            ;;
        --no-py)
            install_py=false
            ;;
        *)
            nprocessors=$arg
#            echo "nproc: $arg $nprocessors"
            ;;
    esac
done

##### End command line arg parsing #####
#echo "nproc: $nprocessors"

echo " Installing GGG"
echo " Using GGGPATH =" $GGGPATH
if [ ! $GGGPATH/install == `pwd` ] ; then
   echo " Your current directory: `pwd`"
   echo " Does not match the GGGPATH install directory."
   read -p " Continue? (Y/N) " req
   if [[ $req == 'Y' || $req == 'y' ]] || $always_yes  ; then
      echo " Continuing..."
   else
      echo " Quitting install. Please change your GGGPATH."
      exit 1
   fi
fi

$GGGPATH/linelist/download_linelists.py

if [[ $? != 0 ]]; then
  echo "An error occurred downloading the linelists (see the previous lines)."
  echo "Correct this error and rerun master.sh. Quitting install."
  exit
fi

ok=`echo $nprocessors | grep -q "^[0-9]*$" && echo "OK" || echo "Not OK"`
if [ "$ok" != "OK" ] ; then
  nprocessors=1
#  echo " Selecting $nprocessors processors."
fi
echo " You've selected $nprocessors processors."

if $install_py; then
    chmod u+x pymaster.sh
    ./pymaster.sh $pyargs
else
    echo "Not cloning/installing python components"
fi

if [ $? != 0 ]; then
    echo "ERROR: Problem cloning netcdf_writer. ABORTING master.sh"
    exit 1
fi

./compile_ggg.sh
echo " ********************************** "
./i2s_master.sh
echo " ********************************** "

echo " Creating data_part.lst if it doesn't already exist"
if [ ! -e ../config/data_part.lst ]; then 
    echo '    The file data_part.lst does not exist. Creating the file...'
    echo "$GGGPATH/spectra/" > ../config/data_part.lst 
else
    echo '    The file data_part.lst already exists. Doing nothing...'
fi
if [ ! -e ../config/data_part_list_maker.lst ]; then 
    echo "$GGGPATH/spectra/" > ../config/data_part_list_maker.lst 
fi
# If the .men files don't exist, copy the version in the install directory
echo " Copying menus..."
if [ ! -e ../levels/levels.men ]; then 
    cp levels.men ../levels/levels.men
fi
# gnd
if [ ! -e ../runlogs/gnd/runlogs.men ]; then 
    cp gnd_runlogs.men ../runlogs/gnd/runlogs.men 
fi
if [ ! -e ../vmrs/gnd/vmrs.men ]; then 
    cp gnd_vmrs.men ../vmrs/gnd/vmrs.men 
fi
if [ ! -e ../models/gnd/models.men ]; then 
    cp gnd_models.men ../models/gnd/models.men
fi
if [ ! -e ../windows/gnd/windows.men ]; then 
    cp gnd_windows.men ../windows/gnd/windows.men
fi
if [ ! -d ../lse ]; then
    mkdir ../lse
    mkdir ../lse/gnd
    mkdir ../lse/lab
elif [ ! -d ../lse/gnd ]; then
    mkdir ../lse/gnd
    if [ ! -d ../lse/lab ]; then
       mkdir ../lse/lab
    fi
elif [ ! -d ../lse/lab ]; then
    mkdir ../lse/lab
fi
# lab
if [ ! -e ../runlogs/lab/runlogs.men ]; then 
   if [ ! -d ../runlogs/lab ]; then
       mkdir ../runlogs/lab
   fi
    cp lab_runlogs.men ../runlogs/lab/runlogs.men 
fi
if [ ! -e ../vmrs/lab/vmrs.men ]; then 
   if [ ! -d ../vmrs/lab ]; then
       mkdir ../vmrs/lab
   fi
   cp lab_vmrs.men ../vmrs/lab/vmrs.men 
fi
if [ ! -e ../models/lab/models.men ]; then 
   if [ ! -d ../models/lab ]; then
       mkdir ../models/lab
   fi
    cp empty_models.men ../models/lab/models.men
fi
if [ ! -e ../windows/lab/windows.men ]; then 
   if [ ! -d ../windows/lab ]; then
       mkdir ../windows/lab
   fi
    cp lab_windows.men ../windows/lab/windows.men
fi
# orb
if [ ! -e ../runlogs/orb/runlogs.men ]; then 
   if [ ! -d ../runlogs/orb ]; then
       mkdir ../runlogs/orb
   fi 
    cp empty_runlogs.men ../runlogs/orb/runlogs.men 
fi
if [ ! -e ../vmrs/orb/vmrs.men ]; then 
   if [ ! -d ../vmrs/orb ]; then
       mkdir ../vmrs/orb
   fi 
    cp empty_vmrs.men ../vmrs/orb/vmrs.men 
fi
if [ ! -e ../models/orb/models.men ]; then 
   if [ ! -d ../models/orb ]; then
       mkdir ../models/orb
   fi 
    cp empty_models.men ../models/orb/models.men
fi
if [ ! -e ../windows/orb/windows.men ]; then 
   if [ ! -d ../windows/orb ]; then
       mkdir ../windows/orb
   fi 
    cp empty_windows.men ../windows/orb/windows.men
fi
# bal
if [ ! -e ../runlogs/bal/runlogs.men ]; then 
   if [ ! -d ../runlogs/bal ]; then
       mkdir ../runlogs/bal
   fi 
    cp empty_runlogs.men ../runlogs/bal/runlogs.men 
fi
if [ ! -e ../vmrs/bal/vmrs.men ]; then 
   if [ ! -d ../vmrs/bal ]; then
       mkdir ../vmrs/bal
   fi 
    cp empty_vmrs.men ../vmrs/bal/vmrs.men 
fi
if [ ! -e ../models/bal/models.men ]; then 
   if [ ! -d ../models/bal ]; then
       mkdir ../models/bal
   fi 
    cp empty_models.men ../models/bal/models.men
fi
if [ ! -e ../windows/bal/windows.men ]; then 
   if [ ! -d ../windows/bal ]; then
       mkdir ../windows/bal
   fi 
    cp empty_windows.men ../windows/bal/windows.men
fi
# mod_maker.input
if [ ! -e ../src/idl/mod_maker.input ]; then
   cp mod_maker.input ../src/idl/mod_maker.input
   echo "$GGGPATH/ncdf/air.2004.parkfalls.nc" >> ../src/idl/mod_maker.input
   echo "$GGGPATH/ncdf/hgt.2004.parkfalls.nc" >> ../src/idl/mod_maker.input
   echo "$GGGPATH/ncdf/pres.tropp.2004.parkfalls.nc" >> ../src/idl/mod_maker.input
   echo "$GGGPATH/ncdf/shum.2004.parkfalls.nc" >> ../src/idl/mod_maker.input
   echo "" >> ../src/idl/mod_maker.input
   echo "This is an example for mod_maker.input." >> ../src/idl/mod_maker.input
   echo "Change your site abbreviation, lat/lon and filenames as appropriate." >> ../src/idl/mod_maker.input
fi
echo " Finished copying menus..."
# Check to ensure that the ak directory exists; if not, make it
if [ ! -d ../ak ]; then
    mkdir ../ak
fi 
# Check to ensure that the current_results directory exists
if [ ! -d current_results ]; then
    mkdir current_results
fi 
rm -f current_results/*
cd current_results
echo " Running list_maker"
../../bin/list_maker < ../.list_maker.input
mv list_maker.out pa_ggg_benchmark.gnd
echo " Running create_sunrun"
../../bin/create_sunrun < ../.create_sunrun.input
echo " Running create_runlog"
../../bin/create_runlog < ../.create_runlog.input
echo " Running gsetup"
../../bin/gsetup < ../.gsetup.input
echo " "
echo " "
echo "GFIT is fitting some benchmark spectra to"
echo "check that the installation was successful."
echo "This will take ~15 minutes."
time_start=`date +%s`

chmod u+x multiggg.sh
# If you wish to run on N processors, use the line starting with "parallel" and comment out the subsequent line. 
# Be sure to change the number after "-j" to the number of processors you wish to use.
if [ "$nprocessors" -eq 1 ] ; then
 echo "Running multiggg.sh"
 ./multiggg.sh
else
 echo "Running multiggg.sh in parallel"
 parallel -j$nprocessors -t --delay 2 < multiggg.sh
fi

chmod u+x post_processing.sh
./post_processing.sh
cd ../

time_end=`date +%s`
time_exec=`expr $(( ($time_end - $time_start)/60 ))`
time_exec_s=`expr $(( ($time_end - $time_start) - $time_exec*60 ))`
echo " "
echo " Execution time was about $time_exec minutes and $time_exec_s second(s)"
echo " "

./diff_versions current_results benchmark_results
