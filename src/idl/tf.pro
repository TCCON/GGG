 pro tf
 filenm=string(" ")
 path=string(" ")
 for i=1,6 do begin
 print,"Enter spectrum  "
 read,filenm
 print,strlen(filenm)
 print,filenm
 ffile,filenm,ip,partition,path
 print,path
 endfor
 end

