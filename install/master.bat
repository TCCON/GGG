rm ../src/*/*.o
compile_all_source_code.bat
cd work
../../bin/create_runlog < cr.inp
../../bin/gsetup < gs.inp
echo " "
echo "Now running GFIT. This will take ~10 minutes. Please be patient."
multiggg.bat
../../bin/gggavg < av.inp
diff pa_h2o.mav original_results/pa_h2o.mav > differences.out
diff pa_h2o.ray original_results/pa_h2o.ray >> differences.out
diff pa_h2o.vav original_results/pa_h2o.vav >> differences.out
diff pa_h2o.vsw original_results/pa_h2o.vsw >> differences.out
echo " "
ls -lt differences.out
cd ../
echo "Finished."
