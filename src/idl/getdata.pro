pro getdata,nfile,fname,param,header,ncol,nrow,gmissing,buf,val,err
;
; Copies the columns of data for which PARAM(ifile)=HEADER(ifile,kcol)
; from BUF(NELE) to VAL and ERR.
;
;  INPUTS:
;          nfile                    Number of input files
;          fname(nfile)             Names (full path) of input files
;          param(nfile)             Column label to be read from each input file
;          header(nfile,1000)       2-D string array of column labels
;          ncol(nfile)              Number of columns in each file
;          nrow(nfile)              Number of rows in each file
;          gmissing                 Global Value corresponding to missing data
;          buf(nele)                Complete data set (in memory).
;
; OUTPUTS:
;          val(ntot)                Array of values
;          err(ntot)                Array of uncertainties
;
;   NOTES:
; 1)If PARAM(ifile) cannot be found in FNAME(ifile), then the appropriate
;   elements of VAL are filled with GMISSING and ERR with 0.0
;
; 2)If PARAM(ifile)_error cannot be found in FNAME(ifile), then the appropriate
;   elements of ERR are filled with 0.0.
;
;
   krow=long(0)
   kele=long(0)
   for ifile=0,nfile-1 do begin
      pname=strtrim(param(ifile),2)
;     Search for PNAME. If not found, kcol=ncol(ifile)
      for kcol=0,ncol(ifile)-1 do if pname eq header(ifile,kcol) then goto,k1
k1:   if kcol lt ncol(ifile) then begin
;         Search for PNAME_error. If not found, ecol=ncol(ifile)
          for ecol=0,ncol(ifile)-1 do if pname+"_error" eq header(ifile,ecol) then goto,k2
k2:       if ecol lt ncol(ifile) then begin
              for irow=long(1),nrow(ifile) do begin
                 val(krow)=buf(kele+kcol)
                 err(krow)=buf(kele+ecol)
                 krow=krow+1
                 kele=kele+ncol(ifile)
              endfor
          endif else begin
              err(krow:krow+nrow(ifile)-1) = 0.0
              for irow=long(1),nrow(ifile) do begin
                 val(krow)=buf(kele+kcol)
                 krow=krow+1
                 kele=kele+ncol(ifile)
              endfor
          endelse
      endif else begin
          print, format='(a80)',pname+string(' not found in ')+fname(ifile)
          val(krow:krow+nrow(ifile)-1) = gmissing
          err(krow:krow+nrow(ifile)-1) = 0.0
          krow=krow+nrow(ifile)
          kele=kele+nrow(ifile)*ncol(ifile)
      endelse
   endfor  ; ifile=0,nfile-1
end
