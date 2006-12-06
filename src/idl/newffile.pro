pro newffile,filenm,pdat,path,flag

path=string(replicate(32b,30))
openr,unit,'/ggg/mkiv/m4part.lst',/get_lun
if pdat.npart eq 0 then begin
  i=0
  repeat begin
    readf,unit,path
    if strmid(path,0,1) ne ':' then i=i+1
  endrep until eof(unit) ne 0
  pdat.npart=i
  pdat=create_struct(pdat,'partition',strarr(pdat.npart))
  point_lun,unit,0
  i=0
  repeat begin
    readf,unit,path
    if strmid(path,0,1) ne ':' then begin
      pdat.partition[i]=path
      i=i+1
    endif
  endrep until eof(unit) ne 0
  close,unit,/all
  free_lun,unit
endif

flag=0
lpart=0
for lpart=0, pdat.npart-1 do begin
  path=string(pdat.partition[pdat.ipart],filenm)
  len=strlen(findfile(path))
  if len[0] gt 0 then return
  pdat.ipart=((pdat.ipart+1) mod (pdat.npart))
endfor
if lpart eq pdat.npart then pdat.ipart=22 $
  else pdat.ipart=pdat.ipart-1
flag=1
return
end
