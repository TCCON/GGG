echo " Deleting .o files"
rm -f ../src/*/*.o
echo " Deleting old executables"
rm -f ../bin/*
cp .gfit_md5sum ../bin/gfit_md5sum
echo " Compiling all programs"
#make -f ../src/*/Makefile
for dir in `ls -1 ../src/ `; do
    if [ -e ../src/${dir}/Makefile ]; then
      make -f ../src/${dir}/Makefile
    fi
done
echo " Deleting .o files"
rm -f ../src/*/*.o
rename .exe "" ../bin/*.exe
cd current_results
rm -f *
echo " Running create_runlog"
../../bin/create_runlog < ../.create_runlog.input
echo " Running gsetup"
../../bin/gsetup < ../.gsetup.input
echo " "
echo " "
echo "GFIT is running on some benchmark spectra"
echo "to check that the installation was successful."
echo "This will take ~15 minutes."
chmod u+x multiggg.sh
./multiggg.sh
chmod u+x post_processing.sh
./post_processing.sh
cd ../
diff --brief --new-file current_results/pa_ggg_benchmark.mav benchmark_results/pa_ggg_benchmark.mav > differences.out
diff --brief --new-file current_results/pa_ggg_benchmark.ray benchmark_results/pa_ggg_benchmark.ray >> differences.out
diff --brief --new-file current_results/pa_ggg_benchmark.vav benchmark_results/pa_ggg_benchmark.vav >> differences.out
diff --brief --new-file current_results/pa_ggg_benchmark.vsw benchmark_results/pa_ggg_benchmark.vsw >> differences.out
echo " "
ls -lt differences.out
echo "If the size of the 'differences.out' file is zero (line above),"
echo "then the installation has completed successfully."
