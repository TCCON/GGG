pro t
mfile=long(10)
fname = strarr(mfile)
buf=fltarr(1000000)
sss =string('xxxxxxxxxxxxxxxx')
for nfile=0,mfile-1 do begin
    print,format='($,"Enter filename ",i1)',nfile
    read, sss
    if (strlen(sss) le 0 ) then goto, ss
    fname(nfile)=strtrim(sss,2)
    openr, unitr, fname(nfile), /get_lun
    readf, unitr, nlhead, nfmt
    close, unitr
    free_lun, unitr
endfor
print, format='("If you have more files, increase parameter MFILE")'
;
ss:

mcol=1000
header    =  strarr(nfile,mcol)
ncol      =  lonarr(nfile)
nrow      =  lonarr(nfile)
ntot=long(0)
nele=long(0)
gmissing=long(-999)
read2mem,nfile,fname,header,ncol,nrow,gmissing,buf
print, '                            Filename            NCOL     NROW     NELE'
for ifile=0,nfile-1 do begin
   ntot=ntot+nrow(ifile)
   nele=nele+nrow(ifile)*ncol(ifile)
   print,format='(a42,3i9)',fname(ifile),ncol(ifile),nrow(ifile),nrow(ifile)*ncol(ifile)
endfor
print,ntot,nele
buf=buf(0:nele-1)
end
