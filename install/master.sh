echo " Deleting .o files"
rm ../src/*/*.o
echo " Deleting old executables"
rm ../bin/*
echo " Compiling all programs"
#make -f ../src/*/Makefile  ! doesn't work so....
for dir in `ls -1 ../src/ `; do
    if [ -d ../src/${dir} ]; then
        make -f ../src/${dir}/Makefile
    fi
done

rename .exe "" ../bin/*.exe
cd current_results
rm *
echo " Running create_runlog"
../../bin/create_runlog < ../input_files/create_runlog.input
echo " Running gsetup"
../../bin/gsetup < ../input_files/gsetup.input
echo " "
echo " "
echo "GFIT is running."
echo "This will take ~15 minutes."
chmod +x multiggg.sh
multiggg.sh
../../bin/collate_results < ../input_files/collate_results.input
../../bin/average_results < ../input_files/average_results.input
../../bin/apply_airmass_correction < ../input_files/apply_airmass_correction.input
../../bin/apply_insitu_correction < ../input_files/apply_insitu_correction.input
cd ../
diff current_results/pa_h2o.mav original_results/pa_h2o.mav > differences.out
diff current_results/pa_h2o.ray original_results/pa_h2o.ray >> differences.out
diff current_results/pa_h2o.vav original_results/pa_h2o.vav >> differences.out
diff current_results/pa_h2o.vsw original_results/pa_h2o.vsw >> differences.out
echo " "
ls -lt differences.out
echo "If the size of the 'differences.out' file is zero (line above),"
echo "then the installation has completed successfully."
