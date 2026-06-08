#!/bin/bash
echo " Testing I2S. This should take a minute or less to complete."
today=`date +"%F %H:%M:%S"`
todayf=`date +"%Y%m%d%H%M%S"`
time_start=`date +%s`
#echo $today > compile_messages.i2s.out
echo $today > $GGGPATH/install/differences.i2s.out
#echo " Compiling i2s."
#cd $GGGPATH/src/i2s/
#echo " Compiling i2s:" >> $GGGPATH/install/compile_messages.i2s.out
#make clean >> $GGGPATH/install/compile_messages.i2s.out 2>&1
#make >> $GGGPATH/install/compile_messages.i2s.out 2>&1
#rm i2s.o
#cd $GGGPATH/src/spec_diff/
#echo " Compiling spec_diff:" >> $GGGPATH/install/compile_messages.i2s.out
#make clean >> $GGGPATH/install/compile_messages.i2s.out 2>&1
#make >> $GGGPATH/install/compile_messages.i2s.out 2>&1
#rm spec_diff.o
cd $GGGPATH/src/i2s/
if [ -e $GGGPATH/src/i2s/spectra/wg20090206saebaa.001 ] ; then
   echo " Deleting Wollongong example spectra and their headers."
   rm $GGGPATH/src/i2s/spectra/wg20090206saeba?.00[1-8]
   if [ -e $GGGPATH/src/i2s/spectra/wg20090206saebaa.001.hdr ] ; then
      rm $GGGPATH/src/i2s/spectra/wg20090206saeba?.00[1-8].hdr
   fi
fi
if [ -e $GGGPATH/src/i2s/spectra/pa20041222saaaaa.001 ] ; then
   echo " Deleting Park Falls example spectra and their headers."
   rm $GGGPATH/src/i2s/spectra/pa20041222saaaa?.00[12]
   if [ -e $GGGPATH/src/i2s/spectra/pa20041222saaaaa.001.hdr ] ; then
      rm $GGGPATH/src/i2s/spectra/pa20041222saaaa?.00[12].hdr
   fi
fi
#echo " Failed to compile:" `grep -E 'Stop|Error' $GGGPATH/install/compile_messages.i2s.out | wc -l` "program(s)."
echo " Running i2s using the Wollongong example input file."
cd $GGGPATH/src/i2s/
$GGGPATH/bin/i2s opus-i2s.example.in > opus-i2s.example.$todayf.i2s.out 2>&1
echo " Running i2s using the Park Falls example input file."
$GGGPATH/bin/i2s slice-i2s.example.in > slice-i2s.example.$todayf.i2s.out 2>&1
#cd $GGGPATH/src/i2s/spectra/
echo " Evaluating differences in the files and headers of the Wollongong and Park Falls i2s spectra." >> $GGGPATH/install/differences.i2s.out
if [ ! -e spectra/wg20090206saebaa.001 ] ; then
  echo '   The Wollongong i2s test failed to write spectra.'
  echo '   The Wollongong i2s test failed to write spectra.' >> $GGGPATH/install/differences.i2s.out
fi
if [ ! -e spectra/pa20041222saaaaa.001 ] ; then
  echo '   The slice-i2s test failed to write spectra.'
  echo '   The slice-i2s test failed to write spectra.' >> $GGGPATH/install/differences.i2s.out
fi

echo " Spectrum differences:" >> $GGGPATH/install/differences.i2s.out
$GGGPATH/src/i2s/spec_diff.sh >> $GGGPATH/install/differences.i2s.out
echo " Spectrum header differences:" >> $GGGPATH/install/differences.i2s.out
$GGGPATH/src/i2s/head_diff.sh >> $GGGPATH/install/differences.i2s.out

# Copy the newly-created spectra into the $GGGPATH/spectra directory for processing by gfit.
# The files are renamed to match the numbering they would have had if the entire days were
# processed. This also renames the spectra to match the GGG2014 runlog.
echo " Copying spectra."
cp $GGGPATH/src/i2s/spectra/pa20040721saaaaa.001 $GGGPATH/spectra/pa20040721saaaaa.043
cp $GGGPATH/src/i2s/spectra/pa20040721saaaab.001 $GGGPATH/spectra/pa20040721saaaab.043
cp $GGGPATH/src/i2s/spectra/pa20040721saaaaa.003 $GGGPATH/spectra/pa20040721saaaaa.119
cp $GGGPATH/src/i2s/spectra/pa20040721saaaab.003 $GGGPATH/spectra/pa20040721saaaab.119
cp $GGGPATH/src/i2s/spectra/pa20041222saaaaa.001 $GGGPATH/spectra/pa20041222saaaaa.019
cp $GGGPATH/src/i2s/spectra/pa20041222saaaab.001 $GGGPATH/spectra/pa20041222saaaab.019
cp $GGGPATH/src/i2s/spectra/pa20041222saaaaa.002 $GGGPATH/spectra/pa20041222saaaaa.020
cp $GGGPATH/src/i2s/spectra/pa20041222saaaab.002 $GGGPATH/spectra/pa20041222saaaab.020

cd $GGGPATH/install/

echo " All done. The results have been written to differences.i2s.out."
time_end=`date +%s`
time_exec=`expr $(( ($time_end - $time_start)/60 ))`
time_exec_s=`expr $(( ($time_end - $time_start) - $time_exec*60 ))`
echo " Execution time was about $time_exec minute(s) and $time_exec_s second(s)"
