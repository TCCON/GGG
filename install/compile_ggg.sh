echo " Deleting .o files"
rm -f ../src/*/*.o
echo " Deleting old executables"
rm -f ../bin/*
cp .gfit_md5sum ../bin/gfit_md5sum
echo " Compiling all programs"
#make -f ../src/*/Makefile
date > compile_messages.out
for dir in `ls -1 ../src/ `; do
    if [ -e ../src/${dir}/Makefile ]; then
      echo "     "${dir}
      make -f ../src/${dir}/Makefile >> compile_messages.out 2>&1
    fi
done
echo " Failed to compile:" `grep -E 'Stop|Error' compile_messages.out | wc -l` "program(s)."
echo " Deleting .o files"
rm -f ../src/*/*.o
echo " Looking for .exe files to rename"
if [ -e ../bin/*.exe ]; then
  echo "     Removing .exe extensions"
  rename .exe "" ../bin/*.exe
else
  echo "     No .exe files found"
fi
