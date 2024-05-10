#! /bin/csh

# We should be in $GGGPATH/src/i2s/
foreach i ( `echo "pa20041222saaaaa.001 pa20041222saaaab.001 pa20041222saaaaa.002 pa20041222saaaab.002 wg20090206saebaa.001 wg20090206saebab.001 wg20090206saebaa.002 wg20090206saebab.002"` )
  echo $i
  ../../utils/OpusHdr/OpusHdr spectra/$i > spectra/$i.hdr
  awk '{print $1}' spectra/benchmark/$i.hdr > tlabnch
  awk '{print $1}' spectra/$i.hdr > tlanew
  grep -v --file=tlabnch spectra/$i.hdr > newtlas
  grep -v --file=tlanew  spectra/benchmark/$i.hdr > oldtlas
  if ( -z newtlas ) then
  cat spectra/$i.hdr > nonewtlas
  else
  grep -v --file=newtlas spectra/$i.hdr > nonewtlas
  endif
  echo 'Param: Benchmark             New_Value     %Difference'
  bash -c 'sdiff -s <( sort spectra/benchmark/'$i'.hdr) <(sort nonewtlas)' | awk '{if ( $2!=0 && $3=="|" ) {print $1,$2,$5,100.* ( $5-$2 ) /$2"%"} else if ($2!=0 && $4=="|" ) {print $1,$2,$6,100.* ( $6-$2 ) /$2"%"} else {print $1,$2,$5}}'
  awk '{print $1,"       NA         ",$2}'  newtlas 
  awk '{print $1,$2,"       NA         "}'  oldtlas
  rm tlabnch tlanew newtlas oldtlas nonewtlas 
end
