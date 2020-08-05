#!/bin/bash
if [ -z $1 ]; then
  nproc=1
else
  nproc=$1
fi
cp .gfit_md5sum ../bin/gfit_md5sum
echo " Compiling all programs and tidying up .o files and old executables"
date > compile_messages.out
time_start=`date +%s`
for dir in `ls -1 ../src/ `; do
    if [ -e ../src/${dir}/Makefile ]; then
      echo "     "${dir} | tee -a compile_messages.out
      make clean -f ../src/${dir}/Makefile >> compile_messages.out 2>&1
      make -j$nproc -f ../src/${dir}/Makefile >> compile_messages.out 2>&1
    fi
done
#there are subroutines found in create_sunrun
for dir in `ls -1 ../src/create_sunrun/ `; do
    if [ -e ../src/create_sunrun/${dir}/Makefile ]; then
      echo "     "create_sunrun/${dir} | tee -a compile_messages.out
      make clean -f ../src/create_sunrun/${dir}/Makefile >> compile_messages.out 2>&1
      make -j$nproc -f ../src/create_sunrun/${dir}/Makefile >> compile_messages.out 2>&1
    fi
done
echo " Failed to properly compile:" `grep -E 'Stop|Error' compile_messages.out | wc -l` "program(s)."
echo " Deleting .o files"
rm -f ../src/*/*.o
rm -f ../src/create_sunrun/*/*.o
echo " Looking for .exe files to rename"
if [ -e ../bin/*.exe ]; then
  echo "     Removing .exe extensions"
  rename .exe "" ../bin/*.exe
else
  echo "     No .exe files found"
fi

time_end=`date +%s`
time_exec=`expr $(( ($time_end - $time_start)/60 ))`
time_exec_s=`expr $(( ($time_end - $time_start) - $time_exec*60 ))`
echo " Compilation time: $time_exec_s s"
