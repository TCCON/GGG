echo " Deleting .o files"
/bin/rm ../src/*/*.o
echo " Deleting old executables"
/bin/rm ../bin/*
echo " Compiling all programs"
make -f  ../src/gfit/Makefile
make -f  ../src/avg_ker/Makefile
make -f  ../src/spec_avg/Makefile
make -f  ../src/list_maker/Makefile
make -f  ../src/dayav/Makefile
make -f  ../src/collate_results/Makefile
make -f  ../src/average_results/Makefile
make -f  ../src/create_runlog/Makefile
make -f  ../src/summarize_linelist/Makefile
make -f  ../src/vll/Makefile
make -f  ../src/gsetup/Makefile
make -f  ../src/extract_pth/Makefile
make -f  ../src/ss2ames/Makefile
make -f  ../src/bin2asc/Makefile
cd current_results
/bin/rm *
echo " Running create_runlog"
../../bin/create_runlog < ../create_runlog.input
echo " Running gsetup"
../../bin/gsetup < ../gsetup.input
echo " "
echo " "
echo "GFIT is running."
echo "This will take ~15 minutes."
echo "Time for a tea break."
chmod +x multiggg.bat
multiggg.bat
../../bin/collate_results < ../collate_results.input
../../bin/average_results < ../average_results.input
cd ../
diff current_results/pa_h2o.mav original_results/pa_h2o.mav > differences.out
diff current_results/pa_h2o.ray original_results/pa_h2o.ray >> differences.out
diff current_results/pa_h2o.vav original_results/pa_h2o.vav >> differences.out
diff current_results/pa_h2o.vsw original_results/pa_h2o.vsw >> differences.out
echo " "
ls -lt differences.out
echo "If the size of the 'differences.out' file is zero (line above),"
echo "then the installation has completed successfully."
