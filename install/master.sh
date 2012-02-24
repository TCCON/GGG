echo " Installing GGG"
echo " Using GGGPATH =" $GGGPATH

./compile_ggg.sh

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
# Check to ensure that the current_results directory exists
if [ ! -d current_results ]; then
    mkdir current_results
fi 
rm -f current_results/*
cd current_results
echo " Running create_runlog"
../../bin/create_runlog < ../.create_runlog.input
echo " Running gsetup"
../../bin/gsetup < ../.gsetup.input
echo " "
echo " "
echo "GFIT is fitting some benchmark spectra to"
echo "check that the installation was successful."
echo "This will take ~15 minutes."
chmod u+x multiggg.sh
./multiggg.sh
chmod u+x post_processing.sh
./post_processing.sh
cd ../

echo " Computing and displaying any differences from the benchmark installation"
diff --brief --new-file current_results/pa_ggg_benchmark.mav benchmark_results/pa_ggg_benchmark.mav > differences.out
diff --brief --new-file current_results/pa_ggg_benchmark.ray benchmark_results/pa_ggg_benchmark.ray >> differences.out
diff --brief --new-file current_results/pa_ggg_benchmark.vav benchmark_results/pa_ggg_benchmark.vav >> differences.out
diff --brief --new-file current_results/pa_ggg_benchmark.vsw benchmark_results/pa_ggg_benchmark.vsw >> differences.out
echo " "
ls -lt differences.out
echo "If the size of the 'differences.out' file is zero (line above),"
echo "then the installation has completed successfully."

echo "Check difference between current_results/pa_ggg_benchmark.mav and benchmark_results/pa_ggg_benchmark.mav" >> differences.out
./ggg_mav_diff.pl current_results/pa_ggg_benchmark.mav benchmark_results/pa_ggg_benchmark.mav >> differences.out
echo "Check difference between current_results/pa_ggg_benchmark.ray and benchmark_results/pa_ggg_benchmark.ray" >> differences.out
./gggdiff.pl current_results/pa_ggg_benchmark.ray benchmark_results/pa_ggg_benchmark.ray >> differences.out
echo "Check difference between current_results/pa_ggg_benchmark.vav and benchmark_results/pa_ggg_benchmark.vav" >> differences.out
./gggdiff.pl current_results/pa_ggg_benchmark.vav benchmark_results/pa_ggg_benchmark.vav >> differences.out
echo "Check difference between current_results/pa_ggg_benchmark.vsw and benchmark_results/pa_ggg_benchmark.vsw" >> differences.out
./gggdiff.pl current_results/pa_ggg_benchmark.vsw benchmark_results/pa_ggg_benchmark.vsw >> differences.out
