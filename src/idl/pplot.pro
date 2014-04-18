pro pplot,datarray,ntgas,npoints,nbox,rmax,targ,xunit,xtxt,ytxt,text,ntxt
;
;  INPUTS:
;       datarray(ntgas+3,npoints)
;       ntgas                      number of target gases
;       npoints                    number of spectral points
;       nbox                       1=no residuals; 2=plot residuals
;       targ                       Index of target gas to be plotted
;       xunit                      =0 cm-1   1=nm
;       xtxt                       x-co-ordinate of text
;       ytxt                       y-co-ordinate of text
;       text                       text to be written
;       ntxt                       Number of text strings to be written
;
;  OUTPUTS:   none
;
;  Implementation Notes:
; 1)  datarray(0,*)       =  Spectral Frequencies
;     datarray(1,*)       =  Measured Spectrum
;     datarray(2,*)       =  Calculated Total Spectrum
;     datarray(3,*)       =  Calculated Spectrum of first target gas
;     datarray(2+i,*)     =  Calculated Spectrum of i'th target gas
;     datarray(2+ntgas,*) =  Calculated Spectrum of last target gas
;
; 2)  text(i) i=0,ntgas-1 =  Target Gas Labels
;     text(ntgas)         =  Title
;     text(ntgas+1)       =  Corner Text (e.g. pfit_000.ps  24-Jul-98/11:04:32)
;     text(ntgas+2+i)     =  User-entered Text Caption (i=0,1,2,3...)
;
; 3)  targ = -1       Only Measured (M) and Calculated (C) spectra plotted
;     targ =  0       0'th target gas transmittance plotted along with M and C
;     targ =  1       1'st target gas transmittance plotted along with M and C
;     targ =  2       2'nd target gas transmittance plotted along with M and C
;     targ =  3       3'rd target gas transmittance plotted along with M and C
;     targ = ntgas-1  Non-target gas transmittance plotted along with M and C
;     targ = ntgas    All gas transmittances plotted along with M and C
;
freq=datarray(0,0:npoints-1)
meas=datarray(1,0:npoints-1)
calc=datarray(2,0:npoints-1)
residual=(meas-calc)
wt=fltarr(npoints)
wt=1/(0.000025+calc*calc/4)  ;  assumes error of 0.5% + calc/2
print,format='(a21,a14,f7.3,a3,f7.3)',strmid(text(ntgas),0,21),": Delta ZOFF =",total(residual*wt)/total(wt)," +-",1/sqrt(total(wt))


!p.charsize=2.5
!p.charsize=2
!p.thick=2
!p.thick=3
!x.thick=2
!y.thick=2
!p.charthick=2
xboxmin=0.14
xboxmax=0.93
yboxmin=0.10
yboxmax=0.90
reswid=0.16 ; Width of the residual box


if nbox eq 1 then begin
   !p.multi=[0,0,1,0,0]
endif else begin
   !p.multi=[0,0,2,0,0]
endelse
!p.position=[xboxmin,yboxmin,xboxmax,yboxmax-reswid*(nbox-1)]
!x.style=1
!x.tickname=''
!x.ticklen=0.02
if xunit eq 0 then begin
!xtitle='!6Frequency (cm!u-1!n)'
endif else begin
!xtitle='!6Wavelength (nm)'
endelse 
!y.style=1
!y.minor=0
!ytitle='!6Transmittance'
plot,freq,calc,linestyle=0
;plot,freq,calc/(datarray(1+ntgas,0:npoints-1)*datarray(3,0:npoints-1)),linestyle=0
oplot,freq,meas,psym=4,symsize=0.5
!p.psym=0
if targ ge ntgas then begin
   tg2=ntgas-1
   tg1=0
endif else begin
   tg2=targ
   tg1=max([tg2,0])
endelse
print,text
xyouts,0.51,0.92,text(ntgas),/normal,alignment=0.5
for i=tg2,tg1,-1 do begin
;   if(i eq 0) then begin
;     kcolor=190
;     kcolor=200
;   endif else if (i eq 1) then begin
;     kcolor=135
;     kcolor=165
;   endif else begin
;     kcolor=190-135*i/(ntgas-1) ; i=0 red; i=ntgas-1 blue
;   endelse
;   kcolor=250-180*i/(ntgas-1) ; i=0 red; i=ntgas-1 blue
;   kcolor=250*(1.0-0.75*i/(ntgas-1)) ; i=0 red; i=ntgas-1 blue
   kcolor=250*(1.0-0.9*i/ntgas) ; i=0 red; i=ntgas-1 blue
;   kcolor=3*kcolor/4
;   if kcolor gt 140 then kcolor = kcolor+30 
   print, i,tg1,tg2,kcolor,ntgas
;   oplot,freq(0:npoints-1),datarray(3+i,0:npoints-1)*datarray(2+ntgas,0:npoints-1),linestyle=0,color=kcolor
   oplot,freq(0:npoints-1),datarray(3+i,0:npoints-1),linestyle=0,color=kcolor
   xyouts,xtxt(i),ytxt(i),text(i),color=kcolor
endfor
;
; Plot calculated spectrum without 1'st target gas contribution
;oplot,freq(0:npoints-1),calc/datarray(3+tg1,0:npoints-1),linestyle=3
;
; Write filename and time in lower left corner
; xyouts,0,0,text(ntgas+1),alignment=0.0,charsize=1.,/normal
;
; Write text captions
for itxt=ntgas+2,ntxt-1 do  xyouts,xtxt(itxt),ytxt(itxt),text(itxt)
;
if nbox ge 2 then begin
   !p.position=[xboxmin,yboxmax-reswid+0.01,xboxmax,yboxmax]
;   !p.position=[0.09,0.76,0.93,0.90]
   !xtitle=''
   !x.tickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' ']
   !x.ticklen=0.10
   !y.minor=-1
;   !ytitle='!6% Residual'
   !ytitle=''
   !y.range=[-rmax,rmax]
   plot,freq,100*residual,linestyle=0
   oplot,freq,(freq*0),linestyle=0
endif
end

