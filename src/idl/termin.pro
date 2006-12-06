pro termin,ttp,termtp,color,idl_device
; initialize variables not passed in as arguments
; (The arguments take on the type of those passed by the calling program)
tchar=''
; get the terminal environment variable
trmnam=strtrim(getenv('TERM'))
; if a null if returned, the TERM variable was unknown.  If this happens
; exit the routine with proper return codes.
if trmnam eq '' then begin
  ttp=-2
  idl_device=''
  return
endif
get_lun,unit
;openr,unit,getenv('IDL_DIR')+'/examples/mark4/common/termid'
openr,unit,getenv('IDL_DIR')+'/examples/common/termid'
readf,unit,nlines
for i=1,nlines do begin
  readf,unit,tchar,termtp,ttp,color,idl_device,$
  format='(a8,1x,a4,2x,i4,2x,a7,2x,a3)'
  tchar=strtrim(tchar)
  if tchar eq trmnam then goto,term_found
  skip_test:
  end
  ttp=-1
  idl_device=''
  term_found:
  close,unit
  idl_device=strtrim(idl_device)
end
