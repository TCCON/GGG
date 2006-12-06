pro slf,yobs,nobs,y0,grad
;
;  Fits a straight line of the form   Yobs(i) = y0 + grad*i
;  Inputs:
;     YOBS(NOBS)  array of data points
;     NOBS        number of data points
;
;  Outputs:
;     Y0          offset of fitted line (with respect to 0'th point)
;     GRAD        gradient of fitted line

tt=0.0
ti=0.0
tii=0.0
tiy=0.0
ty=0.0
for i=long(0),long(nobs-1) do begin
   tt=tt+1
   ti=ti+i
   tii=tii+i*i
   ty=ty+yobs(i)
   tiy=tiy+i*yobs(i)
endfor
denom=ti*ti-tt*tii
y0=(tiy*ti-ty*tii)/denom
grad=(ti*ty-tt*tiy)/denom
end
