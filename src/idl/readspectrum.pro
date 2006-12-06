pro readspectrum,glenfiddich,nus,nue,opd,ifirst,ilast,graw,possp,bpw,$
apo,intrp,freq,yobs,nobs,iflag

iflag=0
res=0.0
rect=0.0
ipmax=6010
print,'readspectrunm'

; print,format='("spectrum= ",a30)',glenfiddich
; print,format='("frqcen= ",f8.3," wid= ",f6.3)',frqcen,wid
; print,format='("nue= ",f9.4," nus= ",f9.4)',nus,nue
; print,format='("ifirst= ",i8," ilast= ",i8)',ifirst,ilast
; print,format='("opd= ",f6.2," possp= ",i8," bpw= ",i3)',opd,possp,bpw
; print,format='("raw grid spacing= ",e14.7)',graw

resn=0.5/opd
if res(0) gt resn(0) then resn=res
resnog=resn/graw
rectog=rect/graw
nk=fix(18*resnog)
if ((apo(0) eq 0) and (intrp(0) eq 1)) then nk=0
ns=fix(2*nk+1)
 print,format='("nk= ",i8," ns= ",i8)',nk,ns

m1=long(nus/graw)+1
m2=long(nue/graw)
 print,format='("m1= ",i8," m2= ",i8)',m1,m2
if graw lt 0. then begin
  xx=m1
  m1=m2
  m2=xx
endif
if m1(0) gt ilast(0) or m2(0) lt ifirst(0) then begin
  iflag=1
  print,format='((a30)," only encompasses ",(f9.3)," to ",(f9.3)," cm-1")',$
  glenfiddich(0),ifirst*graw,ilast*graw
  return
endif
if m1(0)+nk(0) gt ilast(0) or m2(0)-nk(0) lt ifirst(0) then begin
  iflag=1
  print,format='("interval is too close to edge of file to allow convolution")'
 return
endif
if m1(0)-nk(0) lt ifirst(0) then m1=nk+ifirst
if m2(0)+nk(0) gt ilast(0) then m2=ilast-nk
 print,format='("m1= ",i8," m2= ",i8)',m1,m2

openr,unit,glenfiddich,/get_lun
if apo(0) gt 0 or intrp(0) gt 1 then begin
  nsize=m2-m1+ns
  if nsize(0) le 0 then begin
    close,unit
    free_lun,unit
    iflag=1
    return
  endif
  if nsize(0) lt ns(0) then begin
    close,unit
    free_lun,unit
    iflag=1
    return
  endif
  iskip=possp+abs(bpw)*(m1-nk-ifirst)
 print,format='("nsize:",i8," offset:",i8)',nsize,iskip
  case abs(bpw) of
     2: begin
          image=assoc(unit,intarr(nsize(0)),iskip(0))
          spect=image(0)
          if bpw lt 0 then byteorder,spect,/sswap
          ; if ((strpos(glenfiddich,'phg',0) gt 0) or $
          ;     (strpos(glenfiddich,'pin',0) gt 0)) then $
          ;     spect=spect/15000.
        end
     4: begin
          image=assoc(unit,fltarr(nsize(0)),iskip(0))
          spect=image(0)
          if bpw lt 0 then byteorder,spect,/lswap
        end
  endcase
  i11=long(nus/(graw/intrp))+1
  if i11(0) le intrp(0)*(m1(0)-1) then i11=intrp*m1
  i22=long(nue/(graw/intrp))
  if i22(0) ge intrp(0)*(m2(0)+1) then i22=intrp*m2
  print,format='("i11= ",i8," i22= ",i8)',i11,i22
  nobs=i22-i11+1
  print,format='("nobs= ",i8)',nobs
  yobs=fltarr(nobs(0))
  freq=dblarr(nobs(0))
  delwav=graw/intrp
  ; print,format='("delwav= ",e13.7)',delwav
  for q=long(0),long((nobs(0)-1)) do begin
    freq(q)=delwav*(q+i11)
  endfor
  nustrt=i11*delwav
  off=0.0
  a=dblarr(ipmax)
  profwl,apo,intrp*nk,intrp*resnog,intrp*rectog,intrp*off,float(intrp),a,ipmax
  test=fix(intrp)
  test2=intrp*(ns-1)
  for m=0, test-1 do begin
    a(test2+m)=0.
  endfor
  k=long(0)
  test=fix(intrp)
  jj=long(test*m1-i11+1)
  ; print,format='("jj= ",i8)',jj
  for ij=long(0),long(nobs(0)-1) do begin
    dp=0.
    for ix=long(0),long(ns-2) do begin
      dp=dp+a(jj-1+ix*test)*spect(k+ix)
    endfor
    yobs(ij)=dp
    if jj(0) le 1 then begin
      jj=test
      k=k+1
    endif else begin
      jj=jj-1
    endelse
  endfor
  ; print,yobs
endif else begin
  nsize=m2+1-m1
  if nsize(0) le 0 then begin
    close,unit
    free_lun,unit
    iflag=1
    return
  endif
  iskip=abs(bpw)*(m1-ifirst)+possp
  ifirstpoint=ifirst+(1./abs(bpw))*(iskip-possp)-1
  ilastpoint=ifirst+(1./abs(bpw))*(iskip-possp)+(nsize-2)
  ; print,format='("nodp: nsize:",i8," offset:",i8," abs(bpw)= ",i4)',$
  ; nsize,iskip,abs(bpw)
  ; print,format='("start:",i8,"  end:",i8)',ifirstpoint,ilastpoint
  case abs(bpw) of
    2: begin
         image=assoc(unit,intarr(nsize(0)),iskip(0))
         yobs=image(0)
         if bpw lt 0 then byteorder,yobs,/sswap
         ; if ((strpos(glenfiddich,'phg',0) gt 0) or $
         ;     (strpos(glenfiddich,'pin',0) gt 0)) then $
         ;     yobs=yobs/15000.
       end
    4: begin
         image=assoc(unit,fltarr(nsize(0)),iskip(0))
         yobs=image(0)
         if bpw lt 0 then byteorder,yobs,/lswap
       end
  endcase
  freq=dblarr(nsize(0))
  for q=long(0),long(nsize-1) do begin
    freq(q)=graw*(q+m1)
  endfor
endelse
close,unit
free_lun,unit
end
