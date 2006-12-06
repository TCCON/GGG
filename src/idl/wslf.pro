pro wslf,nrow,vx,ex,vy,ey,a,ea,b,eb,chi2on
;  Fits a straight line of the form
;                       y = A.x + B                                  (1)
;  to the data in vectors vx and vy, which have uncertainties ex and ey
;  respectively. Note that the above definition of A and B is the opposite
;  to that used in Numberical recipes (page 504).
;
;  This is achieved by minimizing the quantity
;     CHI2 = SUM { [vy(i)-A.vx(i)-B]**2/[ey(i)**2+(A*ex(i))**2] }    (2)
;  Note that this is much more complicated than the usual case in which
;  ex is assumed to be zero (which makes the equation linear in A by
;  removing it from the denominator). This non-linearity in A necessitates
;  iteration.
;
; INPUTS:
;    NROW      I*4  The number of data points
;    VX(NROW)  R*4  X-Values
;    EX(NROW)  R*4  X-Errors
;    VY(NROW)  R*4  Y-Values
;    EY(NROW)  R*4  Y-Errors
;
; OUTPUTS:
;    A         R*4  Gradient of fitted line
;    EA        R*4  Uncertainty in GRAD
;    B         R*4  Y-Intercept of fitted line
;    EB        R*4  Uncertainty in OFF
;    CHI2ON    R*4  Average factor by which deviations of the data points
;                   from the fitted line exceed their uncertainties.
;
;
; First transform variables so that they all lie in the range -1 to +1
vx_=dblarr(nrow)
ex_=dblarr(nrow)
vy_=dblarr(nrow)
ey_=dblarr(nrow)
xdif=(max(vx(0:nrow-1))-min(vx(0:nrow-1)))/2
xbar=(max(vx(0:nrow-1))+min(vx(0:nrow-1)))/2
ydif=(max(vy(0:nrow-1))-min(vy(0:nrow-1)))/2
ybar=(max(vy(0:nrow-1))+min(vy(0:nrow-1)))/2
   vx_=(vx-xbar)/xdif
   ex_=ex/xdif
   vy_=(vy-ybar)/ydif
   ey_=ey/ydif
tiny=1.d-70
mit=long(19)
a = 0.0
b = 0.0
iflag=long(0)
i=long(0)
iter=long(0)
while (iter lt mit  and  iflag lt 4 ) do begin
   tw  = 0.0d0
   tu  = 0.0
   twx = 0.0
   tux = 0.0
   twy = 0.0
   tuy = 0.0
   for i=long(0),nrow-1 do begin
      ey2 = (1.0d0*ey(i))^2
      ex2 = (1.0d0*ex(i))^2
      dd=ey2+a*ex2*a
      if(dd le tiny) then begin
        w = tiny
        u = tiny*vx(i)
      endif else begin
        w = 1./dd
        u = w*(vx(i)*ey2+a*(vy(i)-b)*ex2)*w
      endelse
      tw  = tw + w
      tu  = tu + u
      twx = twx + w*vx(i)
      tux = tux + u*vx(i)
      twy = twy + w*vy(i)
      tuy = tuy + u*vy(i)
   endfor
   denom=tw*tux-tu*twx
   if (denom le 0.0 ) then  print, ' WSLF:  denom =< 0 ',tw,tux,tu,twx
   awas=a
   bwas=b
   a=(tw*tuy-tu*twy)/denom
   b=(tux*twy-twx*tuy)/denom
   if(abs(a-awas) lt  0.000001*(1.+abs(a))) then iflag=iflag+1
   if(abs(b-bwas) lt  0.000001*(1.+abs(b))) then iflag=iflag+1
   iter=iter+1
endwhile
if(iflag lt 4) then  print,' WSLF: Failed to converge'
ea=sqrt(tw/denom)
eb=sqrt(tux/denom)
;
; Convert back to orignal units.
;ea=ea*ydif/xdif
;eb=sqrt((ydif*eb)^2+(xbar*ea)^2)
;b=ybar+ydif*(b-a*xbar/xdif)
;a=a*ydif/xdif
;
;  Compute Chi-Square
tot=0.0d0
for i=long(0),nrow-1  do begin
   ey2 = ey(i)
   ey2 = ey2^2
   ex2 = ex(i)
   ex2 = ex2^2
   dd=ey2+ex2*a*a
   if(dd le tiny) then begin
       w = tiny
   endif else begin
       w = 1./dd
   endelse
   res = vy(i)-a*vx(i)-b
   tot = tot + w*res^2
endfor
chi2on=sqrt(tot/(nrow-2))
return
end
