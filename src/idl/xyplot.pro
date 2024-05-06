pro xyplot
; General purpose 2-D plotting program.
;
; Version History
version = ' amesplot  V1.3.1   6-Feb-98' ; increased mcol to 1000
version = ' amesplot  V1.3.2  10-Mar-98' ; renames input files to proj.nam
version = ' amesplot  V1.3.3  17-Mar-98' ; does multi-panel plots
version = ' amesplot  V2.0.0  30-Mar-98' ; first reads all data to memory
version = ' amesplot  V2.0.1  30-Apr-98' ; increased max # input files to 40
version = ' amesplot  V2.0.2  11-May-98' ; added 'txt' command.
version = ' amesplot  V2.0.3  17-May-98' ; added 'captsize' keyword.
version = ' amesplot  V2.0.4  28-May-98' ; merged "rf" & "pr" commands.
version = ' amesplot  V2.0.5   8-Jun-98' ; Added NTIMES to SCX & SCY commands
version = ' amesplot  V2.0.6  17-Jun-98' ; included cornertext as text(0)
version = ' amesplot  V2.0.7  17-Jun-98' ; No longer uses SUBSTR for TXT or CP.
version = ' amesplot  V2.0.8  20-Aug-98' ; Introduced "rc?" and "help" commands
version = ' amesplot  V2.0.9  11-Sep-98' ; Added "LABELS:" keyword.
version = ' amesplot  V2.1.0  11-Sep-98' ; Added x2b,b2x,y2b,b2y,dyb,dxb commands
version = '   xyplot  V3.0.0   4-Feb-99' ; Introduced use of registers  
version = '   xyplot  V3.1.0   8-Jul-99' ; Abbreviated vd#, vm# and va# commands
version = '   xyplot  V3.1.1  26-Jul-99' ; Cleaned up the code.
version = '   xyplot  V3.2.0  20-Sep-99' ; New "ce" command (Copy Error).
version = '   xyplot  V3.2.1  21-Sep-99' ; Fixed bug in pythag.
version = '   xyplot  V3.3.0  28-Sep-99' ; New "np" command (New Page).
version = '   xyplot  V3.3.1  21-Oct-99' ; Added "ln" cmd (Natural Logarithm).
version = '   xyplot  V3.3.2   3-Jan-00' ; Fixed minor bugs in sd command.
version = '   xyplot  V3.3.3  16-Jan-00' ; Added 'lf' (straight line fit).
version = '   xyplot  V3.3.4   6-Mar-00' ; Added 'ex' (exponent).
version = '   xyplot  V3.3.5  29-Mar-00' ; Added 'vs' (vector subtract).
version = '   xyplot  V3.3.6   7-Apr-00' ; Removed /Times flag from device call
version = '   xyplot  V3.3.7  19-Apr-00' ; psn=1 so first plot is now _001.ps
version = '   xyplot  V3.3.8  25-Oct-00' ; Fixed bug in Ames format 2160
version = '   xyplot  V3.3.9  12-Dec-00' ; Introduced "ticklen" keyword to .sym
version = '   xyplot  V3.4.0  23-Jan-01' ; Allowed NAUXV=0 in Ames 2010 format.
version = '   xyplot  V3.4.1  23-May-02' ; Can now skip NDSC header line in Ames files
version = '   xyplot  V3.5.0  11-Nov-03' ;  Made cc commanding more efficient
version = '   xyplot  V3.5.1  11-Sep-06' ;  
version = '   xyplot  V3.5.2  11-Sep-06' ;  
version = '   xyplot  V3.5.3  21-Dec-07' ;  
version = '   xyplot  V3.54  2018-04-17' ;  


; The following code should already be inside your ~/.idl/startup.pro
; In case not, it is included below.
if(!version.os_family eq 'unix') then device, true_color=24
window, /free, /pixmap, colors=-10
wdelete, !d.window
device, retain=2, decomposed=0, set_character_size=[10,12]
device, get_visual_depth=depth
print, 'Display depth = ', depth
print, 'Color table size = ', !d.table_size

; Initialize stuff
isp=0
ipanel=0
flen=1
ntxt=1
ntimes=1.0
psn=1
kslf=0

sval=string(' ')
anno=string(' ')
sss =string(' ')
ccc =string(' ')
outputfile=string(' ')
proj=string('xyplot')

print
print,version
print
mfile=50  ; Maximum number of input files
fname = strarr(mfile)

;openw, units, 'av.out', /get_lun
openw, unitw, 'xyplot.log', /get_lun
printf,unitw,'xyplot'

readfiles:
nccreg=0
; Prompt user for names of input data files.
for nfile=0,mfile-1 do begin
    ljj=0
    while ljj(0) le 0 do begin
    print,format='($,"Enter name of input data file (CR if done) ",i2)',nfile
    read, sss
    printf, unitw, sss
    if (strlen(sss) le 0 ) then goto, ss
    ljj=strlen(findfile(sss))
    if ljj(0) le 0 then print, ' Unable to find file. Please try again.'
    endwhile
;
;  Check that the file can be opened and the first header line can be read.
    fname(nfile)=strtrim(sss,2)
    openr, unitr, fname(nfile), /get_lun
    readf, unitr, ccc
    if strmid(ccc,0,4) eq 'TOON' then readf, unitr, ccc    ; skip NDSC header
    reads, ccc, nlhead, nfmt
    close, unitr
    free_lun, unitr
endfor
print, format='("If you have more files, increase parameter MFILE")'

ss:
if nfile le 0 then begin
   print, 'No input data files were been found'
   goto, readfiles
endif

print,'Current Project Name: '+proj
print,'Enter Project Name (CR to leave unchanged)'
read, sss
printf,unitw,sss
if strlen(sss) gt 0 then   proj=strtrim(sss,2)
ljj=strlen(findfile(proj+'.sym'))
if ljj(0) gt 0 then begin
   symfile=proj+'.sym'
endif else begin
   symfile=string('xyplot.sym')
   ljj=strlen(findfile(symfile))
   if ljj(0) le 0 then $
      symfile=string(getenv('GGGPATH')+'/src/idl/xyplot.sym')
   print,' File not found: '+proj+'.sym'
endelse
print, symfile

param     =  strarr(nfile)
ncol      =  lonarr(nfile)
nrow      =  lonarr(nfile)

print, ' Reading entire contents of data files into memory.....'
buf = fltarr(160000000)
mcol=1200
header    =  strarr(nfile,mcol)
gmissing   =  -999.
gmissing   =  9.8765E+35

read2mem,nfile,fname,gmissing,header,ncol,nrow,buf

ntot=long(0)
nele=long(0)
print, '                            Filename           NCOL     NROW     NELE'
for ifile=0,nfile-1 do begin
   ntot=ntot+nrow(ifile)
   nele=nele+nrow(ifile)*ncol(ifile)
   print,format='(a42,3i9)',fname(ifile),ncol(ifile),nrow(ifile),nrow(ifile)*ncol(ifile)
endfor
buf=buf(0:nele-1)

xx2plot = fltarr(ntot)
ex2plot = fltarr(ntot)
yy2plot = fltarr(ntot)
ey2plot = fltarr(ntot)
aa2plot = fltarr(ntot)
bb2plot = fltarr(ntot)

msym=400  ;  maximum number of different symbols
color_index  = lonarr(msym)
caption= strarr(msym)
xcp    = fltarr(msym)
ycp    = fltarr(msym)

mtxt=20  ;  maximum number of permitted text strings.
xtxt   = fltarr(mtxt)
ytxt   = fltarr(mtxt)
text   = strarr(mtxt)
text(0)='x'

;  mreg is number of registers available for storage of plotable data
;  i.e. data already fetched from the memory file dump (BUF) by GETDATA.
;  The first two registers (0 & 1) are the data that are plotted (x & y)
;  The upper registers are used for temporary storage and for controlling
;  modifications to the plot.
regname=string('xyabcdefg')
mreg  = strlen(regname)
ccmin  =  fltarr(mreg,nfile)
ccmax  =  fltarr(mreg,nfile)
ccindex=  lonarr(mreg)
iccindex=  lonarr(mreg)
ccindex(*)=-1
iccindex(*)=-1
vmin  = fltarr(mreg)
vmax  = fltarr(mreg)
cr    = intarr(mreg)
title = strarr(nfile,mreg+1)
value = fltarr(ntot,mreg+1)
error = fltarr(ntot,mreg+1)
result = fltarr(ntot)

;value(*,0)=findgen(ntot)

repeat begin
   print
   banner=string('Register                 File0                    File1                    File2')
   print,strmid(banner,0,10+25*nfile)
   for ireg=0,mreg-1 do  print,format='(2x,a1,2x,3a25)', strmid(regname,ireg,1),title(0:min([2,nfile-1]),ireg(0))
   print,format='($,"Enter rd,sm,sd,sa,vm,va,lb,au,pd,ca,wf,qu ")'
   read, sss
   printf, unitw, sss
   sss=strtrim(sss,2)
   lens=strlen(sss)
   if lens ge 2 then cmd=strmid(sss,0,2)
;  Last character of command string (for lens > 2) is destination vector.
;
;  Extract the names of the destination (kd) and auxiliary (ka) registers.
   kd=mreg
   for kw=2,lens-1 do begin
      ka=kd
      for kk=0,mreg-1 do begin
         if strmid(sss,kw,1) eq strmid(regname,kk,1) then kd=kk
      endfor
   endfor
;
   case cmd of

     'he' : begin  ; Help
	print,'  Recognized XYPLOT Commands include:'
	print,'     rd#  -  reads selected quantity into the #-register'
	print,'     rc#  -  takes the reciprocal of the contents of the #-register'
	print,'     lg#  -  takes base 10 logarithm of contents of the #-register'
	print,'     ln#  -  takes natural logarithm of contents of the #-register'
	print,'     ex#  -  takes exponential of contents of the #-register'
	print,'     sr#  -  takes square-root of contents of the #-register'
	print,'     co#  -  takes the cosine of the contents of the #-register'
	print,'     sm#  -  multiply the #-register by a user-supplied scalar'
	print,'     sd#  -  divide the #-register by a user-supplied scalar'
	print,'     sa#  -  add a user-supplied scalar to the #-register'
        print,'     ss#  -  Selective Scalar Add to the #-register'
	print,'     va&# -  adds the &-register to the #-register'
	print,'     vs&# -  subtracts the &-register from the #-register'
	print,'     vm&# -  multiplies the #-register by the &-register'
	print,'     vd&# -  divides the #-register by the &-register'
	print,'     cp&# -  copies the &-register value and error to the #-register'
	print,'     ce&# -  copies the &-register error to the #-register value'
	print,'     cr#  -  manually Change Range of #-axis'
	print,'     au#  -  automatically select range of #-axis'
	print,'     cc#  -  only plot points for which min<#-register<max
	print,'     lb#  -  label for #-axis'
	print,'     wf   -  writes output file of x- and y-registers'
	print,'     rf   -  returns to the read file(s) point'
	print,'     rr   -  re-reads files'
	print,'     np   -  forces New Page to allow incomplete multipanel plots'
	print,'     ca   -  writes a caption onto plot'
	print,'     tx   -  writes text onto plot'
	print,'     pd   -  makes a screen plot'
	print,'     hc   -  makes a hard copy'
	print,'     qu   -  quits program'
	print,' where # denotes a destination (i.e. input/output) register'
	print,' and & denotes an auxilliary (i.e. input) register,'
	print,' Both # and & represent x, y, a, b, c, d, e, f'
        end

     'lg' : begin  ; Logarithm (base 10)
          bad=where(value(*,kd) le 0.0  or value(*,kd) eq gmissing, nmiss)
          error(*,kd) = error(*,kd) / value(*,kd) / 2.3026
          value(*,kd) = alog10( value(*,kd) )
          if(nmiss gt 0) then value(bad,kd)=gmissing
          title(*,kd)="log!d10!n["+title(*,kd)+"]"
          end

     'ln' : begin  ; Natural Logarithm (base e)
          bad=where(value(*,kd) le 0.0  or value(*,kd) eq gmissing, nmiss)
          error(*,kd) = error(*,kd) / value(*,kd)
          value(*,kd) = alog( value(*,kd) )
          if(nmiss gt 0) then value(bad,kd)=gmissing
          title(*,kd)="ln["+title(*,kd)+"]"
          end

     'ex' : begin  ; Exponentiate (the inverse of 'ln')
          bad=where(value(*,kd) gt 88.0  or  value(*,kd) eq gmissing, nmiss)
          value(*,kd) = exp( value(*,kd) )
          error(*,kd) = error(*,kd) * value(*,kd)
          if(nmiss gt 0) then value(bad,kd)=gmissing
          title(*,kd)="exp["+title(*,kd)+"]"
          end

     'sr' : begin  ;  Square Root
          bad=where(value(*,kd) lt 0.0  or  value(*,kd) eq gmissing, nmiss)
          value(*,kd) = sqrt( value(*,kd) )
          error(*,kd) = error(*,kd) / value(*,kd) / 2
          if(nmiss gt 0) then value(bad,kd)=gmissing
          title(*,kd)="sqrt["+title(*,kd)+"]"
          end

     'si' : begin  ; Sine
          value(*,kd) = value(*,kd)*3.14159265/180.
          error(*,kd) = error(*,kd)*3.14159265/180.
          value(*,kd) = sin( value(*,kd) )
          error(*,kd) = error(*,kd)*cos(value(*,1))
          title(*,kd)="sine["+title(*,kd)+"]"
          end

     'co' : begin  ; Cosine
          value(*,kd) = value(*,kd)*3.14159265/180.
          error(*,kd) = error(*,kd)*3.14159265/180.
          value(*,kd) = cos( value(*,kd) )
          error(*,kd) = error(*,kd)*sin(value(*,1))
          title(*,kd)="cosine["+title(*,kd)+"]"
          end

     'va' : begin  ; Vector Add

          if ka ge mreg then begin  ;  no auxiliary register specified
          getparam, unitw, header, ncol, fname, param, nfile
          title(*,ka) = param
          cp2reg,nfile,fname,param,header,ncol,nrow,gmissing,buf,ka,value,error
          endif

          print,format='($," Number of times to be added")'
          read,sss
          printf,unitw,sss
          reads,sss,ntimes
          if ntimes gt 1 then begin
             title(*,kd) = "("+title(*,kd) + "+" + sss+"*"+title(*,ka)+")"
          endif else begin
             title(*,kd) = "("+title(*,kd) + "+" + title(*,ka)+")"
          endelse
          bad=where(value(*,ka) eq gmissing or value(*,kd) eq gmissing, nmiss)
          value(*,kd) = value(*,kd) + ntimes*value(*,ka)
          pythag, error(*,kd), ntimes*error(*,ka), result, ntot
          error(*,kd)=result
          if(nmiss gt 0) then value(bad,kd)=gmissing
          end

     'vs' : begin  ; Vector Subtract

          if ka ge mreg then begin  ;  no auxiliary register specified
            getparam, unitw, header, ncol, fname, param, nfile
            title(*,ka) = param
          cp2reg,nfile,fname,param,header,ncol,nrow,gmissing,buf,ka,value,error
          endif
          title(*,kd) = "("+title(*,kd) + "-" + title(*,ka)+")"
          bad=where(value(*,ka) eq gmissing or value(*,kd) eq gmissing, nmiss)
          value(*,kd) = value(*,kd) - value(*,ka)
          pythag, error(*,kd), error(*,ka), result, ntot
          error(*,kd)=result
          if(nmiss gt 0) then value(bad,kd)=gmissing
          end

     'sm' : begin  ; Scalar Multiply
          khi=-1
          bad=where(value(*,kd) eq gmissing, nmiss)
          for ifile=0,nfile-1 do begin
             klo=khi+1
             khi=khi+nrow(ifile)
             print,format='($," Scalar multiplier to be applied to ",a)',fname(ifile)+">"+title(ifile,kd)
             read, sval
             printf,unitw,sval
             reads, sval, xmul
             value(klo:khi,kd) = value(klo:khi,kd) * xmul
             error(klo:khi,kd) = error(klo:khi,kd) * xmul
             title(ifile,kd)=sval+"*["+title(ifile,kd)+"]"
          endfor
          if (nmiss gt 0) then value(bad,kd) = gmissing
          end                   

     'sd' : begin  ; Scalar Divide
          khi=-1
          bad=where(value(*,kd) eq gmissing, nmiss)
          for ifile=0,nfile-1 do begin
             klo=khi+1
             khi=khi+nrow(ifile)
             print,format='($," Scalar by which to divide ",a)',fname(ifile)+">"+title(ifile,kd)
             read, sval
             printf,unitw,sval
             reads, sval, xdiv
             value(klo:khi,kd) = value(klo:khi,kd) / xdiv
             error(klo:khi,kd) = error(klo:khi,kd) / xdiv
             title(ifile,kd)="["+title(ifile,kd)+"]/"+sval
          endfor
          if (nmiss gt 0) then value(bad,kd) = gmissing
          end                   

     'sa' : begin  ; Scalar Add
          khi=-1
          bad=where(value(*,kd) eq gmissing, nmiss)
          for ifile=0,nfile-1 do begin
             klo=khi+1
             khi=khi+nrow(ifile)
             print,format='($,"Enter scalar to be added to ",a)',fname(ifile)+">"+title(ifile,kd)
             read, sval
             printf,unitw,sval
             reads, sval, xadd
             if(nrow(ifile) gt 0) then value(klo:khi,kd) = value(klo:khi,kd) + xadd
             title(ifile,kd)=sval+"+"+title(ifile,kd)
          endfor
          if (nmiss gt 0) then value(bad,kd) = gmissing
          end

     'ra' : begin  ; Ramp (missing x-values)
          khi=-1
          for ifile=0,nfile-1 do begin
             xx=0l
             klo=khi+1
             khi=khi+nrow(ifile)
             for jj=klo,khi do begin
                value(jj:jj,kd) = xx
                xx=xx+1
             endfor
          title(ifile,kd)="Ramp"
          endfor
          end

     'lr' : begin  ; Limit Range
          khi=-1
          bad=where(value(*,kd) eq gmissing, nmiss)
          print,format='($,"Enter Min/Max limits",a)'
          read, sval
          printf,unitw,sval
          reads, sval, zmin,zmax
          print,zmin,zmax
          value(where(value(*,kd) lt zmin),kd) = zmin
          value(where(value(*,kd) gt zmax),kd) = zmax
          if (nmiss gt 0) then value(bad,kd) = gmissing
          end

     'ss' : begin  ; Selective Scalar Add
          khi=-1
          bad=where(value(*,kd) eq gmissing, nmiss)
          print,format='($,"Enter scalar to be added and threshold",a)'
          read, sval
          printf,unitw,sval
          reads, sval, xadd, thresh
          add=where(value(*,kd) lt thresh, icount)
          print,icount
          if (icount gt 0) then value(add,kd)=value(add,kd)+xadd
          if (nmiss gt 0) then value(bad,kd) = gmissing
          end

     'lb' : begin  ; Label
          print,format='($,"Enter annotation:")'
          read, anno
          printf,unitw,anno
          title(*,kd) = anno
          end

     'np' : begin  ;  Force New Page
          device, /close
          ipanel=0
          psn=psn+1
          end

     'sp' : begin  ;  skip panel
          isp=isp+1
          ipanel=ipanel+1
          psn=psn+1
          end

     'au' : cr(kd) = 0  ; automatic range control

     'rr' : begin  ;  Re-read the input data file(s)
          read2mem,nfile,fname,gmissing,header,ncol,nrow,buf
          for kk=0,mreg-1 do begin
            param = title(*,kk) 
            np=strlen(param)
            if (np(0) gt 0) then cp2reg,nfile,fname,param,header,ncol,nrow,gmissing,buf,kk,value,error
          endfor
          end

     'rf' : goto,readfiles ; Read new files

     'rd' : begin    ; Read into a register
          if kd lt mreg then begin
             getparam, unitw, header, ncol, fname, param, nfile
             cp2reg,nfile,fname,param,header,ncol,nrow,gmissing,buf,kd,value,error
             title(*,kd) = param
          endif else begin
             print, sss+' ??   Unknown register'
          endelse
          end

     'cc' : begin  ; Constrain data plotted
          if kd lt mreg then begin
             for ifile=0,nfile-1 do begin
                print,format='($,"Enter plot range (min,max) for ",a)',fname(ifile)+">"+title(ifile,kd)
                read,zmin,zmax
                printf,unitw,zmin,zmax
                ccmin(kd,ifile)=zmin
                ccmax(kd,ifile)=zmax
             endfor
             if ccindex(kd) lt 0  then begin ; It's a new constraint
                ccindex(kd)=nccreg
                iccindex(nccreg)=kd
                nccreg=nccreg+1
             endif
          endif else begin
             print,'invalid destination register'
          endelse
          end

     'in' : begin  ;  Integer Conversion 
          title(*,kd) = "INT("+title(*,kd)+")"
          value(*,kd)=long(value(*,kd))
          end

     'ab' : begin  ;  Absolute Value
          title(*,kd) = "|"+title(*,kd)+"|"
          value(*,kd)=abs(value(*,kd))
          end


     'vm' : begin  ; Vector Multiply

          if ka ge mreg then begin  ;  no auxiliary register specified
             getparam, unitw, header, ncol, fname, param, nfile
             title(*,ka) = param
             cp2reg,nfile,fname,param,header,ncol,nrow,gmissing,buf,ka,value,error
          endif
;
          title(*,kd) = title(*,kd)+"*"+title(*,ka)
          bad=where(value(*,ka) eq gmissing or value(*,kd) eq gmissing, nmiss)
          pythag, error(*,kd)*value(*,ka), value(*,kd)*error(*,ka), result, ntot
          error(*,kd)=result
          value(*,kd) = value(*,kd)*value(*,ka)
          if(nmiss gt 0) then value(bad,kd)=gmissing
          end

     'dp' : begin  ; Dot product of two vectors
          bad=where(value(*,ka) eq gmissing or value(*,kd) eq gmissing, nmiss)
          if nmiss eq 0 then begin
            khi=-1
            for ifile=0,nfile-1 do begin
              klo=khi+1
              khi=khi+nrow(ifile)
              dotprod=0.0
              dperr2=0.0
              print,klo,khi
              for kk=klo,khi  do begin
                dotprod=dotprod+value(kk,kd)*value(kk,ka)
                dperr2 =dperr2+(error(kk,kd)*value(kk,ka))^2+(value(kk,kd)*error(kk,ka))^2
              endfor
              print, ifile,dotprod,sqrt(dperr2)
            endfor
          endif
          end
;

     'rc' : begin  ; takes the reciprocal of the vector
          title(*,kd) = '1 / '+title(*,kd)
          bad=where(value(*,kd) eq gmissing or value(*,kd) eq 0, nmiss)
          error(*,kd) = error(*,kd)/value(*,kd)/value(*,kd)
          value(*,kd) = 1/value(*,kd)
          if nmiss gt 0  then  value(bad,kd) = gmissing
          end

     'vp' : begin  ;  computes vapor pressure (of H2O)
          title(*,kd) = 'SVP_H2O'
          bad=where(value(*,kd) eq gmissing or value(*,kd) eq 0, nmiss)
          error(*,kd) = error(*,kd)*(6111.72784/value(*,kd)+0.15215)/value(*,kd)
          value(*,kd) = exp(23.33086-6111.72784/value(*,kd)+0.15215*alog(value(*,kd)))
          error(*,kd) = error(*,kd)*value(*,kd)
          if nmiss gt 0  then  value(bad,kd) = gmissing
          end

     'vd' : begin  ; Vector Divide
          if ka ge mreg then begin  ;  no auxiliary register specified
             getparam, unitw, header, ncol, fname, param, nfile
             title(*,ka) = param
             cp2reg,nfile,fname,param,header,ncol,nrow,gmissing,buf,ka,value,error
          endif

          title(*,kd) = title(*,kd)+'/'+title(*,ka)
          bad=where(value(*,ka) eq gmissing or value(*,ka) eq 0 or value(*,kd) eq gmissing, nmiss)
          if nmiss gt 0  then  value(bad,ka) = gmissing
          value(*,kd) = value(*,kd)/value(*,ka)
;          print,error(*,kd)/value(*,ka), value(*,kd)*error(*,ka)/value(*,ka)
          pythag, error(*,kd)/value(*,ka), value(*,kd)*error(*,ka)/value(*,ka), result, ntot
          error(*,kd)=result
          if nmiss gt 0  then  value(bad,kd) = gmissing
          end

     'cp' : begin ; copy contents and attributes of one register to another register
          if lens eq 4 then begin
             value(*,kd) = value(*,ka)
             error(*,kd) = error(*,ka)
             title(*,kd) = title(*,ka)
             vmin(kd)    = vmin(ka)
             vmax(kd)    = vmax(ka)
             cr(kd)      = cr(ka)
          endif else begin
             print,"Try again: cp command requires 2 arguments"
          endelse
          end

     'ce' : begin ; copy uncertainties of one register to another register
          if lens eq 4 then begin
             value(*,kd) = error(*,ka)
             error(*,kd) = 0
             title(*,kd) = title(*,ka)+'_error'
;             cr(kd)      = 0
          endif else begin
             print,"Try again: ce command requires 2 arguments"
          endelse
          end

     'pd' : begin  ;  Make a screen plot
          !xtitle=title(0,0)
          !ytitle=title(0,1)
          !p.title=title(0,2)
          set_plot, 'x'
          npp=long(0)  ;  total number of plottable points
          nsym=long(0)
          jrow=long(0)
          bbmin=+9.99e+36
          bbmax=-9.99e+36
          for ifile=0,nfile-1 do begin
             temp=-999999
;                   color_index(nsym)=npp
;                   nsym = nsym + 1
             for i = long(1), nrow(ifile) do begin
                if temp ne fix(value(jrow,2)) then begin  ; a-register sets color switches
                   temp = fix(value(jrow,2))
                   color_index(nsym)=npp
                   nsym = nsym + 1
                   totx=0.0
                   toty=0.0
                   totex=0.0
                   totey=0.0
                endif
                if value(jrow,0) eq gmissing then  goto, skip1
                if value(jrow,1) eq gmissing then  goto, skip1
;                for ireg=0,mreg-1 do begin
                for ireg=0,nccreg-1 do begin
                   kc=iccindex(ireg)
                   if value(jrow,kc) lt ccmin(kc,ifile) then goto, skip1
                   if value(jrow,kc) gt ccmax(kc,ifile) then goto, skip1
                endfor ;  ireg=0,mreg-1 
                xx2plot(npp)=value(jrow,0)
                ex2plot(npp)=error(jrow,0)
                yy2plot(npp)=value(jrow,1)
                ey2plot(npp)=error(jrow,1)
                aa2plot(npp)=value(jrow,2)
                bb2plot(npp)=value(jrow,3)
                totx=totx+value(jrow,0)
                totex=totex+error(jrow,0)^2
                toty=toty+value(jrow,1)
                totey=totey+error(jrow,1)^2
                if (bb2plot(npp) lt bbmin) then bbmin=bb2plot(npp)
                if (bb2plot(npp) gt bbmax) then bbmax=bb2plot(npp)
                npp=npp+1
skip1:
                jrow=jrow+1
             endfor  ; i = long(1), nrow(ifile)
             color_index(nsym)=npp
          endfor     ;  ifile=0,nfile-1
          print,'nsym,npp=',nsym,npp

          if( npp le 0) then  begin
             print, " Warning - no plotable points"
          endif else begin

          !x.style = 17  ; force exact x-range
          if cr(0) eq 0 then begin  ; Automatically choose X-Range X
             set_range,[min(xx2plot(0:npp-1)),max(xx2plot(0:npp-1))],xout,xseglen
             !x.range = xout
          endif else begin
             set_range,[vmin(0), vmax(0)],xout,xseglen
             !x.range = [vmin(0), vmax(0)]
          endelse
;          print,'x_range: ',min(xx2plot(0:npp-1)),max(xx2plot(0:npp-1))
;          print,'set_range: xout: ',xout,xseglen

          !y.style = 17  ; force exact y-range
          if cr(1) eq 0 then begin ; Automatically choose X-Range
             set_range,[min(yy2plot(0:npp-1)),max(yy2plot(0:npp-1))],yout,yseglen
             !y.range = yout
          endif else begin
             set_range,[vmin(1), vmax(1)],yout,yseglen
             !y.range = [vmin(1), vmax(1)]
          endelse
;          print,'y_range: ',min(yy2plot(0:npp-1)),max(yy2plot(0:npp-1))
;          print,'set_range: yout: ',yout,yseglen
       endelse
       plotname=proj+'_'+string(format='(i3.3)',psn)+'.ps'
       time=strmid(systime(0),8,2)+'-'+strmid(systime(0),4,3)+'-'+$
       strmid(systime(0),20,4)+'/'+strmid(systime(0),11,8)
;       if(strlen(strtrim(text(0),2)) gt 0) then  text(0)=plotname+'  '+time
       gplot, xx2plot,ex2plot,yy2plot,ey2plot,color_index,xcp,ycp,symfile,caption,nsym,ipanel,npanel,orientation,xtxt,ytxt,text,ntxt
;       hplot,xx2plot,ex2plot,yy2plot,ey2plot,bbmin,bbmax,bb2plot,xcp,ycp,symfile,caption,npp,ipanel,npanel,orientation,xtxt,ytxt,text,ntxt,xseglen,yseglen
       if kslf gt 0 then oplot,xi,yi,psym=0,thick=3  ; Add fitted straight line
       end

     'av' : begin  ; Average  X and Y values
          klo=color_index(0)
          tiny=1.e-16
          for jsym=0,nsym-1 do begin
            khi=color_index(jsym+1)-1
            avx=-999
            avy=-999
            ex=-999
            ey=-999
            if(khi ge klo) then begin
               tx=0.0
               twx=0.0
               ty=0.0
               twy=0.0
               for k=klo,khi do begin
                  if(ex2plot(k) lt tiny) then ex2plot(k)=tiny
                  if(ey2plot(k) lt tiny) then ey2plot(k)=tiny
                  wx=1.0/ex2plot(k)^2
                  wy=1.0/ey2plot(k)^2
                  tx=tx+wx
                  ty=ty+wy
                  twx=twx+wx*xx2plot(k)
                  twy=twy+wy*yy2plot(k)
               endfor
               avx=twx/tx
               avy=twy/ty
               twx=0.0
               twy=0.0
               for k=klo,khi do begin
                  wx=1.0/ex2plot(k)^2
                  wy=1.0/ey2plot(k)^2
                  twx=twx+wx*(xx2plot(k)-avx)^2
                  twy=twy+wy*(yy2plot(k)-avy)^2
               endfor
               ex=sqrt(twx/tx)
               ey=sqrt(twy/ty)
            endif
            print,format='(i5,4f12.3)',khi-klo+1,avx,ex,avy,ey
;            printf,units,format='(i5,4f12.3)',khi-klo+1,avx,ex,avy,ey
            klo=khi+1
          endfor ; jsym=0,nsym-1
          end

     'lf' : begin  ; Straight line fit
          klo=color_index(0)
          for jsym=0,nsym-1 do begin
             khi=color_index(jsym+1)-1
             if(khi ge klo) then begin
                slfbe,khi-klo+1,xx2plot(klo:khi),ex2plot(klo:khi),yy2plot(klo:khi),ey2plot(klo:khi),a,ea,b,eb,chi2on,pcc
                print,format='(a,2i5,5e13.4,f9.5)','lf: slfbe: ',jsym+1,khi-klo+1,a,ea*chi2on,b,eb*chi2on,chi2on,pcc
             endif
             klo=khi+1
          endfor ; jsym=0,nsym-1
          slfbe,klo,xx2plot,ex2plot,yy2plot,ey2plot,a,ea,b,eb,chi2on,pcc
          print,format='(a,i5,5e13.4,f9.5)','lf: slfbe: all: ',klo,a,ea*chi2on,b,eb*chi2on,chi2on,pcc
          line_rectangle_intersect,!x.range,!y.range,a,b,xi,yi,istat
          kslf=1
          end

     'wf' : begin  ; Write ASCII file containing contents of x- and y-registers
          print, ' Enter Name of Output file'
          read, outputfile
          printf, unitw, outputfile
          openw, unitwf, outputfile, /get_lun
          printf, unitwf, '2 5'
          printf, unitwf, title(0,0)+' ',title(0,0)+'_error   ',title(0,1)+' ',title(0,1)+'_error   bb '
          for i=0,npp-1 do begin
          printf, unitwf, format='(5e14.6)',xx2plot(i),ex2plot(i),yy2plot(i),ey2plot(i),bb2plot(i)
          endfor
          close, unitwf
          free_lun, unitwf
          print,' written file '+outputfile
          end

     'ca' : begin  ;  Write symbols and captions onto plot
          for isym=0,nsym-1 do begin
             print,format='($,i1,"/",i1,") Enter Xcap, Ycap, caption")',isym+1,nsym
             read, xd,yd,sss
             printf, unitw, xd,yd,sss
             caption(isym)=strmid(sss,1,999)
             xcp(isym)=xd
             ycp(isym)=yd
          endfor
          end

     'tx' : begin  ;  Write text onto plot.
          for itxt=0,ntxt-1 do print,format='(i2,2e11.4,1x,a)',itxt,xtxt(itxt),ytxt(itxt),text(itxt)
          print,' Enter #N, X, Y, Text'
          read, itxt,xt,yt,sss
          printf, unitw, itxt,xt,yt,sss
          text(itxt)=strmid(sss,1,999)  ; skip  SSS(0) which is the delimiter
          xtxt(itxt)=xt
          ytxt(itxt)=yt
          ntxt=max([itxt+1,ntxt])
          end


     'hc' :  begin  ;  Write current plot to a postscript file
          set_plot, 'ps'
          if ipanel-isp eq 0 then begin  ;  open a new plotfile
;             !p.font=7   ; this is necessary to get the Times font
;             !p.font=1   ; this is necessary to get the TrueType font
             if (orientation eq "portrait") then begin 
                device,/portrait,font_size=10,filename=plotname,/color,xoffset=2.0,yoffset=1.0,xsize=17.4,ysize=23.5
             endif else begin
                device,/landscape,font_size=10,filename=plotname,/color
             endelse
          endif
          gplot,xx2plot,ex2plot,yy2plot,ey2plot,color_index,xcp,ycp,symfile,caption,nsym,ipanel,npanel,orientation,xtxt,ytxt,text,ntxt
;          hplot,xx2plot,ex2plot,yy2plot,ey2plot,bbmin,bbmax,bb2plot,xcp,ycp,symfile,caption,npp,ipanel,npanel,orientation,xtxt,ytxt,text,ntxt,xseglen,yseglen
          if kslf gt 0 then begin
             oplot,xi,yi,psym=0,thick=4  ; plot fitted straight lines
;             oplot,[1.01*xi(0),0.],[yi(0),(yi(0)*xi(1)-xi(0)*yi(1))/(xi(1)-xi(0)),yi(1)],psym=0,thick=4  ; N2O-OCS balloon correlation figure
;             oplot,[0.,0.],[119.,0.],psym=0,thick=4 ; N2O-OCS balloon correlation figure
             text_slf_b= '!9D!6y/!9D!6x ='+string(format='(f5.1)',b)+'!9+!6'+string(format='(f3.1)',4*eb*chi2on)
             text_slf_r= string(format='("PCC =",(f5.2))',pcc)
             xcap=0.975*!x.range(0)+0.025*!x.range(1)
             xyouts,xcap,0.91*!y.range(1),text_slf_b,charsize=1
             xyouts,xcap,0.83*!y.range(1),text_slf_r,charsize=1
             kslf=0
          endif
          ipanel=ipanel+1
          if ntxt gt 1 then  text(1:ntxt-1)=' '
          ntxt=1
          caption(*)=''
          if ipanel eq npanel then begin  ; close plotfile
             device, /close
             ipanel=0
             psn=psn+1
          endif
          end

     'cr' :  begin  ; Set range of plotted values.
          cr:
          on_ioerror,cr
          print,format='($,"Enter min and max range")'
          read,sval
          printf,unitw,sval
          reads,sval,zmin,zmax
          vmin(kd)=zmin
          vmax(kd)=zmax
          cr(kd) = 1
          end

     'fl' :  begin  ; Frame Left
          del=vmax(kd)-vmin(kd)
          vmax(kd)=vmin(kd)
          vmin(kd)=vmax(kd)-del
          end

     'fr' :  begin  ; Frame Right
          del=vmax(kd)-vmin(kd)
          vmin(kd)=vmax(kd)
          vmax(kd)=vmin(kd)+del
          end

     'qu' :  begin  ; Quit.
          if ipanel gt 0 then   device, /close  ; close plotfile
          close, unitw
          free_lun, unitw
          end

     else :  print, 'Command not recognized...'

   endcase

endrep until cmd eq 'qu'
 
end
