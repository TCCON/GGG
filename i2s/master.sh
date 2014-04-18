echo " Installing I2S. This should take a minute or less to complete."
today=`date +"%F %H:%M:%S"`
todayf=`date +"%Y%m%d%H%M%S"`
time_start=`date +%s`
echo $today > compile_messages.out
echo $today > differences.out
echo " Compiling opus-i2s."
cd opus-i2s/
echo " Compiling opus-i2s:" >> ../compile_messages.out
make clean >> ../compile_messages.out 2>&1
make >> ../compile_messages.out 2>&1
rm opus-i2s.o
cd ..
echo " Compiling slice-i2s."
cd slice-i2s/ 
echo "---------------------" >> ../compile_messages.out
echo " Compiling slice-i2s:" >> ../compile_messages.out
make clean >> ../compile_messages.out 2>&1
make >> ../compile_messages.out 2>&1
rm slice-i2s.o
cd ..
echo " Compiling cobs."
cd cobs/ 
echo "---------------------" >> ../compile_messages.out
echo " Compiling cobs:" >> ../compile_messages.out
make clean >> ../compile_messages.out 2>&1
make >> ../compile_messages.out 2>&1
rm cobs.o
cd ..
echo " Cleaning up."
rm comn/*.o
rm opus-comn/*.o
if [ -e opus-i2s/spectra/wg20090206saebaa.001 ] ; then
   echo " Deleting Wollongong (opus-i2s) spectra and their headers."
   rm opus-i2s/spectra/wg20090206saeba?.00[1-8]
   if [ -e opus-i2s/spectra/wg20090206saebaa.001.hdr ] ; then
      rm opus-i2s/spectra/wg20090206saeba?.00[1-8].hdr
   fi
fi
if [ -e slice-i2s/spectra/pa20041222saaaaa.001 ] ; then
   echo " Deleting Park Falls (slice-i2s) spectra and their headers."
   rm slice-i2s/spectra/pa20041222saaaa?.00[12]
   if [ -e slice-i2s/spectra/pa20041222saaaaa.001.hdr ] ; then
      rm slice-i2s/spectra/pa20041222saaaa?.00[12].hdr
   fi
fi
echo " Failed to compile:" `grep -E 'Stop|Error' compile_messages.out | wc -l` "program(s)."
echo " Running opus-i2s using the example input file."
cd opus-i2s/
./opus-i2s opus-i2s.example.in > opus-i2s.example.$todayf.out 2>&1
cd spectra/
echo " Evaluating differences in the files and headers of the opus-i2s spectra." >> ../../differences.out
if [ ! -e wg20090206saebaa.001 ] ; then
  echo '   The opus-i2s test failed to write spectra.'
  echo '   The opus-i2s test failed to write spectra.' >> ../../differences.out
else
  ../../cobs/cobs ../../cobs/opus-i2s.in
  tail -n+4 ../../cobs/opus-i2s.in.out >> ../../differences.out
  for i in benchmark/wg*.0??
    do
      j=`echo $i | awk -F/ '{print $2}'`
      if [ ! -e $j ] ; then
        echo "Error: File $j *does not exist*" >> ../../differences.out
      else
      ../../scripts/OpusHdr $j > $j.hdr
      diff -q $j benchmark/$j > /dev/null
      if [ $? == 0 ] ; then
        echo "Files $j and benchmark/$j are the same" >> ../../differences.out
      else
#        diff --brief --new-file $j benchmark/$j >> ../../differences.out
        diff -q $j.hdr benchmark/$j.hdr > /dev/null
        if [ $? == 0 ] ; then
          echo "Files $j.hdr and benchmark/$j.hdr are the same" >> ../../differences.out
        else
#          diff --brief --new-file $j.hdr benchmark/$j.hdr >> ../../differences.out
          echo "---------" >> ../../differences.out
          echo "$j.hdr benchmark/$j.hdr % diff" >> ../../differences.out
          sdiff -s $j.hdr benchmark/$j.hdr | awk '{if($5!=0){print $1,$2,$5,100.*($2-$5)/$5"%"} else {print $1,$2,$5}}' >> ../../differences.out
        fi
      fi
      fi
  done
fi
cd ../../
echo " Running slice-i2s using the example input file."
cd slice-i2s/
./slice-i2s slice-i2s.example.in > slice-i2s.example.$todayf.out 2>&1
cd spectra/
echo "--------------------------------------------------------------------------" >> ../../differences.out
echo " Evaluating differences in the files and headers of the slice-i2s spectra." >> ../../differences.out
if [ ! -e pa20041222saaaaa.001 ] ; then
  echo '   The slice-i2s test failed to write spectra.'
  echo '   The slice-i2s test failed to write spectra.' >> ../../differences.out
else
  ../../cobs/cobs ../../cobs/slice-i2s.in
  tail -n+4 ../../cobs/slice-i2s.in.out >> ../../differences.out
  for i in benchmark/pa*.0??
    do
      j=`echo $i | awk -F/ '{print $2}'`
      if [ ! -e $j ] ; then
        echo "Error: File $j *does not exist*" >> ../../differences.out
      else
      ../../scripts/OpusHdr $j > $j.hdr
      diff -q $j benchmark/$j > /dev/null
      if [ $? == 0 ] ; then
        echo "Files $j and benchmark/$j are the same" >> ../../differences.out
      else
#        diff --brief --new-file $j benchmark/$j >> ../../differences.out
        diff -q $j.hdr benchmark/$j.hdr > /dev/null
        if [ $? == 0 ] ; then
          echo "Files $j.hdr and benchmark/$j.hdr are the same" >> ../../differences.out
        else
#          diff --brief --new-file $j.hdr benchmark/$j.hdr >> ../../differences.out
          echo "---------" >> ../../differences.out
          echo "$j.hdr benchmark/$j.hdr % diff" >> ../../differences.out
          sdiff -s $j.hdr benchmark/$j.hdr | awk '{if($5!=0){print $1,$2,$5,100.*($2-$5)/$5"%"} else {print $1,$2,$5}}' >> ../../differences.out
        fi
      fi
      fi
  done
fi
cd ../../
echo " All done. The results have been written to differences.out."
time_end=`date +%s`
time_exec=`expr $(( ($time_end - $time_start)/60 ))`
time_exec_s=`expr $(( ($time_end - $time_start) - $time_exec*60 ))`
echo " Execution time was about $time_exec minute(s) and $time_exec_s second(s)"
