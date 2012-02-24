echo " Deleting .o files"
rm -f ../src/*/*.o
echo " Deleting old executables"
rm -f ../bin/*
cp .gfit_md5sum ../bin/gfit_md5sum
echo " Compiling all programs"
#make -f ../src/*/Makefile  ! doesn't work so....
for dir in `ls -1 ../src/ `; do
    if [ -d ../src/${dir} ]; then
        make -f ../src/${dir}/Makefile
    fi
done

rename .exe "" ../bin/*.exe
cd current_results
rm -f *
echo " Running create_runlog"
../../bin/create_runlog < ../.create_runlog.input
echo " Running gsetup"
../../bin/gsetup < ../.gsetup.input
echo " "
echo " "
echo "GFIT is running."
echo "This will take ~15 minutes."
chmod +x multiggg.sh
./multiggg.sh
chmod +x post_processing.sh
./post_processing.sh
cd ../
diff current_results/pa_h2o.mav original_results/pa_h2o.mav > differences.out
diff current_results/pa_h2o.ray original_results/pa_h2o.ray >> differences.out
diff current_results/pa_h2o.vav original_results/pa_h2o.vav >> differences.out
diff current_results/pa_h2o.vsw original_results/pa_h2o.vsw >> differences.out
echo " "
ls -lt differences.out
echo "If the size of the 'differences.out' file is zero (line above),"
echo "then the installation has completed successfully."
