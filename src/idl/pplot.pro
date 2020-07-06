pro pplot,xmin,xmax,ymin,ymax,datarray,ntgas,npoints,nbox,rmax,ktgas,xzo,xunit,xtxt,ytxt,text,ntxt
;
;  INPUTS:
;       datarray(ntgas+4,npoints)
;       ntgas                  number of target gases
;       npoints                number of spectral points
;       nbox                   1=no residuals; 2=plot residuals
;       ktgas                  Index of target gas to be plotted
;       xzo                    zero level offset (transmittance)
;       xunit                  =0 cm-1   1=nm
;       xtxt                   x-co-ordinate of text
;       ytxt                   y-co-ordinate of text
;       text                   text to be written
;       ntxt                   Number of text strings to be written
;
;  OUTPUTS:   none
;
;  Implementation Notes:
; 1)  datarray(0,*)       =  Spectral Frequencies
;     datarray(1,*)       =  Measured Spectral Transmittance
;     datarray(2,*)       =  Calculated Spectral Spectrum
;     datarray(3,*)       =  Calculated Continuum (AU)
;     datarray(4,*)       =  Calculated Transmittance of first target gas
;     datarray(3+i,*)     =  Calculated Transmittance of i'th target gas
;     datarray(3+ntgas,*) =  Calculated Transmittance of last target gas
;
; 2)  text(i) i=0,ntgas-1 =  Target Gas Labels
;     text(ntgas)         =  Title
;     text(ntgas+1)       =  Corner Text (e.g. pfit_000.ps  24-Jul-98/11:04:32)
;     text(ntgas+2+i)     =  User-entered Text Caption (i=0,1,2,3...)
;
; 3)  ktgas = -1       Only Measured (M) and Calculated (C) spectra plotted
;     ktgas =  0       0'th target gas transmittance plotted along with M and C
;     ktgas =  1       1'st target gas transmittance plotted along with M and C
;     ktgas =  2       2'nd target gas transmittance plotted along with M and C
;     ktgas =  3       3'rd target gas transmittance plotted along with M and C
;      .      .        .
;      .      .        .
;     ktgas = ntgas-1  Non-target gas transmittance plotted along with M and C
;     ktgas = ntgas    All gas transmittances plotted along with M and C
;     ktgas = ntgas+1  Cont
;
;  SPTS*TRAN = [Calcul/Cont-XZO]/(1-XZO)
;  Calcul = (SPTS*TRAN*(1-XZO)+XZO)*Cont 

freq=datarray(0,0:npoints-1)
meas=datarray(1,0:npoints-1)
calc=datarray(2,0:npoints-1)
residual=(meas-calc)

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
!xtitle='!6Wavenumber (cm!u-1!n)'
endif else begin
!xtitle='!6Wavelength (nm)'
endelse 
!y.style=1
!y.minor=0
!ytitle='!6Transmittance'

; Define xtickformat
frange=xmax-xmin
ifm=ceil(alog10(max([abs(xmin),abs(xmax)])))
ifr=round(alog10(0.85*abs(frange)+1.E-18))
if ifr gt 0 then xtf='(i'+string(format='(i1)',ifm)+')' else $
xtf='(f'+string(format='(i1)',ifm-ifr+3)+'.'+string(format='(i1)',1-ifr)+')'

print,'xfrange,ifm,ifr,xtf = ',frange,ifm,ifr,xtf

; Define ytickformat
;print,ymin,ymax
yrange=ymax-ymin
ifm=ceil(alog10(max([abs(ymin),abs(ymax)])))
if( ymin lt 0 ) then ifm=ifm+1  ; ' allow space for minus sign
ifr=round(alog10(abs(0.85*yrange)+1e-18))
if ifm gt 9 then begin
   ytf='(i'+string(format='(i2)',ifm)+')'
endif else if ifr gt 0 then begin
   ytf='(i'+string(format='(i1)',ifm)+')' 
endif else begin
    ytf='(f'+string(format='(i1)',ifm-ifr+3)+'.'+string(format='(i1)',1-ifr)+')'
endelse
print,'ymin,ymax,yrange,ifm,ifr,ytf = ',ymin,ymax,yrange,ifm,ifr,ytf


!p.psym=0
if ktgas ge ntgas+1 then begin  ; Continuum
   tg2=-1
   tg1=-1
   !ytitle='!6Signal (A.U.)'
;  Calcul = (SPTS*TRAN*(1-XZO)+XZO)*Cont 
   calc=(calc*(1-xzo)+xzo)*datarray(3,*)
   meas=meas*datarray(3,*)
endif else if ktgas ge ntgas then begin  ; all target gases
   tg2=ntgas-1
   tg1=0
endif else begin  ; one target gas (ktgas)
   tg2=ktgas
   tg1=max([tg2,0])
endelse
plot,freq,calc,linestyle=0,xtickformat=xtf,ytickformat=ytf
oplot,freq,meas,psym=4,symsize=0.5
xyouts,0.51,0.92,text(ntgas),/normal,alignment=0.5
;for i=tg1,tg2 do begin
for i=tg2,tg1,-1 do begin
   kcolor=242*(1.0-0.9*i/ntgas) ; i=0 red; i=ntgas-1 blue
   if i eq 0 then kcolor=248
   if kcolor lt 190 then kcolor=kcolor-15
;   oplot,freq(0:npoints-1),datarray(4+i,0:npoints-1)*datarray(3+ntgas,0:npoints-1),linestyle=0,color=kcolor
   oplot,freq(0:npoints-1),datarray(4+i,0:npoints-1),linestyle=0,color=kcolor
   xyouts,xtxt(i),ytxt(i),text(i),color=kcolor
endfor
;
; Plot calculated spectrum without 1'st target gas contribution
;oplot,freq(0:npoints-1),calc/datarray(4+tg1,0:npoints-1),linestyle=3
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
   !ytitle='!6Residual'
;   !ytitle=''
   !y.range=[-rmax,rmax]
   plot,freq,residual,linestyle=1,psym=4,symsize=0.4
   oplot,freq,residual,linestyle=0
   oplot,freq,(freq*0),linestyle=0
endif
end
