./delete_all_temporary_files.bat
rm ../bin/*
./delete_all_object_files.bat
./compile_all_source_code.bat
../bin/create_sunrun_from_parkfalls_ifs1 < cs.inp
../bin/create_runlog < cr.inp
../bin/gsetup < gs.inp
chmod +x multiggg.bat
./multiggg.bat
../bin/gggavg < av.inp
./check_files_for_differences.bat
