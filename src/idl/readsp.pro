pro readsp,path,nus,nue,opd,ifirst,ilast,graw,possp,bpw,$
apo,intrp,freq,yobs,nout,iflag

iflag=0
ccc = string(' ')
print,'readsp'

; print,format='("spectrum= ",a30)',path
; print,format='("frqcen= ",f8.3," wid= ",f6.3)',frqcen,wid
; print,format='("nue= ",f9.4," nus= ",f9.4)',nus,nue
; print,format='("ifirst= ",i8," ilast= ",i8)',ifirst,ilast
 print,format='("opd= ",f6.2," possp= ",i8," bpw= ",i3)',opd,possp,bpw
; print,format='("raw grid spacing= ",e14.7)',graw

resn=0.5/opd
rect=0.0
resnog=resn/graw
rectog=rect/graw
nhw=fix(18*resnog)
if ((apo(0) eq 0) and (intrp(0) eq 1)) then nhw=0
nop=2*intrp*nhw+1
;print,format='("nhw= ",i8," nop= ",i8)',nhw,nop

m1=long(intrp*nus/graw)+1
m2=long(intrp*nue/graw)
if graw lt 0. then begin
   xx=m1
   m1=m2
   m2=xx
endif
n1=(m1+intrp-1)/intrp-nhw   ; index of first spectral point needed
n2=m2/intrp+nhw             ; index of last spectral point needed
delwav=graw/intrp
if n1(0)+2*nhw gt ilast(0) or n2(0)-2*nhw lt ifirst(0) then begin
  iflag=1
  print,format='((a30)," only encompasses ",(f9.1)," to ",(f9.1)," cm-1")',$
  path(0),ifirst*graw,ilast*graw
  return
endif
if n1(0) lt ifirst(0) then begin
   n1=ifirst
   m1=intrp*(n1+nhw)-intrp+1
endif
if n2(0) gt ilast(0) then begin
   n2=ilast
   m2=intrp*(n2-nhw)
endif
 print, nus,nue,graw
 print,format='("m1= ",i8," m2= ",i8)',m1,m2
nout=m2-m1+1
nin=n2-n1+1
freq=dblarr(nout)
for q=long(0),long((nout(0)-1)) do begin
   freq(q)=delwav*(q+m1)
endfor

ff=0.0
tm=0.0
signal=0.0
openr,unit,path,/get_lun
 iskip=possp+abs(bpw)*(n1-ifirst)
 print,format='("nin:",i8," iskip:",i8)',nin,iskip
 getendian,iend
; iend=-1
  print,'bpw,iend=',bpw,iend(0)
  case abs(bpw) of
     2: begin
          image=assoc(unit,intarr(nin(0)),iskip(0))
          yobs=image(0)
          if bpw*iend(0) lt 0 then byteorder,yobs,/sswap
        end
     4: begin
          image=assoc(unit,fltarr(nin(0)),iskip(0))
          yobs=image(0)
          if bpw*iend(0) lt 0 then byteorder,yobs,/lswap
        end
     7: begin
          readf, unit, nlhead  ;  Skip comment lines
          for k=1,nlhead-1 do begin
             readf, unit, ccc
          endfor
          for i=long(0),long(n1-ifirst) do  begin
             readf, unit,ff,tm,signal
          endfor
          yobs=fltarr(nin(0))
          for i=long(0),long(nin(0))-1 do begin
             readf, unit,ff,tm,signal
             yobs(i)=signal
          endfor
        end
     9: begin
          readf, unit, nlhead  ;  Skip comment lines
          for k=1,nlhead-1 do begin
             readf, unit, ccc
          endfor
          for i=long(0),long(n1-ifirst) do  begin
             readf, unit,ff,signal
          endfor
          yobs=fltarr(nin(0))
          for i=long(0),long(nin(0))-1 do begin
             readf, unit,ff,signal
             yobs(i)=signal
          endfor
        end
  endcase
  for i=long(0),long(nin-1) do begin
   ff=graw*(i+n1)
yobs=float(yobs)
; FTUVS
;  cnt=(1.1+0.15*exp(-((ff-9970)/1070)^2)+3.35*exp(-((ff-11900)/1070)^2))/(0.108+exp(-((ff-8250)/2550)^2))
;  yobs(i)=(0.003*ff/13000+float(yobs(i)))*cnt
;
; Kitt Peak
;  cnt=0.36E-06*(1.00+5.0*exp(-((ff-14650)/800)^2)+2.95*exp(-((ff-12650)/1400)^2)-0.09*exp(-((ff-8850)/500)^2)-0.44*exp(-((ff-10900)/800)^2)+5.5*exp(-((ff-7850)/400)^2)+0.7*exp(-((ff-12050)/1800)^2)+0.2*exp(-((ff-10600)/400)^2)+0.1*exp(-((ff-9400)/150)^2)-0.1*exp(-((ff-11450)/300)^2))
;  yobs(i)=(0.000+float(yobs(i)))*cnt
;  cnt=0.19E-06*(1+0.1*exp(-((ff-7300)/900)^2)+((ff-6400)/1400)^2)
;  yobs(i)=(0.000+float(yobs(i)))*cnt
  endfor

if apo(0) gt 0 or intrp(0) gt 1 then begin
  yout=fltarr(nout)
  off=0.0
  a=fltarr(nop)
  profwl,apo,intrp*nhw,intrp*resnog,intrp*rectog,intrp*off,a
  for i=1,intrp do begin
     dp=0.0
     for j=0,2*nhw-1 do begin
        dp=dp+a(intrp*j+i)
     endfor
     for j=0,2*nhw-1 do begin
        a(intrp*j+i)=a(intrp*j+i)/dp
     endfor
  endfor
  print,nin,nop,nout
  newdec,yobs,nin,a,nop,intrp,yout,nout
  yobs=yout
endif
close,unit
free_lun,unit
end
