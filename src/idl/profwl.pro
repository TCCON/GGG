pro profwl,apo,nk,resnos,rectos,off,a
can=dblarr(1)
t=dblarr(1)
t2=dblarr(1)
t4=dblarr(1)
p=dblarr(1)
q0=dblarr(1)
q1=dblarr(1)
q2=dblarr(1)
q4=dblarr(1)
dd=dblarr(1)
totl=dblarr(1)
c=fltarr(4,4)
c(0,0)=1.0
c(0,1)=0.0
c(0,2)=0.0
c(0,3)=0.0
c(1,0)=0.548000
c(1,1)=-8.33000e-02
c(1,2)=0.535300
c(1,3)=0.0
c(2,0)=0.260000
c(2,1)=-0.154838
c(2,2)=0.894838
c(2,3)=0.0
c(3,0)=9.00000e-02
c(3,1)=0.587500
c(3,2)=0.322500
c(3,3)=0.3225
if (apo(0) gt 3) then apo=3
if (apo(0) lt 0) then apo=0
if (nk le 0) then begin
  a(0)=1.0
  return
endif
;
can(0)=(3.1415926565/resnos)
for ij=0,2*nk(0)-1 do begin
  a(ij)=0.0
  xx=ij-nk(0)
  t(0)=double(abs(can(0)*(xx-off(0))))
  t2(0)=double(t(0)*t(0))
  t4(0)=double(t2(0)*t2(0))
  if (t(0) ge 1.1) then begin
    q0(0)=sin(t(0))/t(0)
    p(0)=cos(t(0))
    q1=3*(q0(0)-p(0))/t2(0)
    q2=-15*((1-3/t2(0))*Q0+3*P/t2(0))/t2(0)
    q4=945*((1-15*(3-7/t2(0))/t2(0))*q0(0)+5*(2-21/t2(0))*p(0)/t2(0))/t4(0)
  endif else begin
    q0(0)=1-t2(0)*(1-t2(0)*(1-t2(0)*(1-t2(0)*(1-t2(0)/110)/72)/42)/20)/6
    q1(0)=1-t2(0)*(1-t2(0)*(1-t2(0)*(1-t2(0)*(1-t2(0)/130)/88)/54)/28)/10
    q2(0)=1-t2(0)*(1-t2(0)*(1-t2(0)*(1-t2(0)*(1-t2(0)/150)/104)/66)/36)/14
    q4(0)=1-t2(0)*(1-t2(0)*(1-t2(0)*(1-t2(0)*(1-t2(0)/170)/136)/90)/52)/22
  endelse
    a(ij)=c(apo,0)*q0(0)+c(apo,1)*q1(0)+c(apo,2)*q2(0)+c(apo,3)*q4(0)
endfor
; weakly apodize slit function
dd(0)=nk(0)*nk(0)
for ij=0,2*nk(0)-1 do begin
  t2(0)=1.0-xx*xx/dd(0)
  a(ij)=a(ij)*t2(0)*t2(0)
endfor
return
end
