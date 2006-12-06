pro twslf
ccc = string(' ')
nlhead = long(0)
ncol = long(0)
mrow = long(100)
nrow = long(0)
vx = fltarr(mrow)
ex = fltarr(mrow)
vy = fltarr(mrow)
ey = fltarr(mrow)
openr, unit, '$GGGPATH/f77/tslf.dat', /get_lun
readf, unit, nlhead, ncol
for i=2,nlhead do  readf, unit, ccc
i=long(0)
while not EOF(unit)  do begin
   readf, unit, av,ae,bv,be
   vx(i)=av
   ex(i)=ae
   vy(i)=bv
   ey(i)=be
   i=i+1
endwhile
close, unit
nrow=i
a=0.0
ea=0.0
b=0.0
eb=0.0
chi2on=0.0
wslf,nrow,vx,ex,vy,ey,a,ea,b,eb,chi2on
print,nrow,a,ea,b,eb,chi2on
end
