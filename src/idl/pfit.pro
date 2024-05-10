pro pfit

; The following code should already be inside your ~/.idl/startup.pro
; In case not, it is included below.
if(!version.os_family eq 'unix') then device, true_color=24
window, /free, /pixmap, colors=-10
wdelete, !d.window
device, retain=2, decomposed=0, set_character_size=[8,10]
device, get_visual_depth=depth
print, 'Display depth = ', depth
print, 'Color table size = ', !d.table_size

ncol=0
;mss=20   ; Maximum number of sub-strings (Freq, Tm  Tc  + 15 target gases)
xunit=1  ; nm
xunit=0  ; cm-1
ymin=0.0
ymax=1.0
color=string('t')
sss=string(' ')
headers=string(' ')
;gases=strarr(mss)
idl_device=string('x')
rlrecord=string('')
occul=string('flt91125.bal')
filenm=string('pfit_001.ps')
filenm=string('')
nbox=2
ktgas=-1
rlim=string("a")
atemp=string("a")
xlolim=string("a")
xuplim=string("a")
ylolim=string("a")
yuplim=string("a")
same_window=0
npoints=long(0)
mspec=long(720000)
sptpath=strarr(mspec)
termtp=string('xwin')
xtxt   = fltarr(80)
ytxt   = fltarr(80)
mtxt   = intarr(80)
text   = strarr(80)
ntxt   = 0
mtxt(*)=0

psn=1
idn=0
ttyp=0
vbar_was=0.
version=string(' PFIT   Version 4.82     2018-04-22    GCT ')
print,version

disk=getenv('GGGPATH')+'/'

;termin,ttyp,termtp,color,idl_device
set_plot,idl_device
window,0,title='pfit.pro',xsize=1100,ysize=512
xmenu,['N-2','N-1 (Previous)','N+0 (Replot)','N+1 (Next)','N+2','Skip','Residual','R_Max','R_Auto',$
  'X_Min','X_Max','X_Auto','X_Left','X_Right','X_Zoom','X_Anti-Zoom',$
  'Y_Min','Y_Max','Y_Auto','Target_Gases','Text','Same_Window','cm-1/nm',$
  'Hardcopy','Quit'],base=b1,button=butn
widget_control,/realize,b1

rtn1:
on_ioerror,rtn1
print,format='($,"enter runlog")'
read,occul
;
;idot=strpos(occul,'.')    ; First instance of '.' in string OCCUL
;
;  Find last instance of '.' in string OCCUL
idot=strlen(occul)
while (strmid(occul,idot,1) ne '.') do begin
   idot=idot-1
endwhile
print,idot+1, strmid(occul,idot+1,1)
;
case strmid(occul,idot+1,1) of
   'g': openr,lunr_rl,string(disk+'runlogs/gnd/'+occul),/get_lun
   'a': openr,lunr_rl,string(disk+'runlogs/air/'+occul),/get_lun
   'b': openr,lunr_rl,string(disk+'runlogs/bal/'+occul),/get_lun
   'o': openr,lunr_rl,string(disk+'runlogs/orb/'+occul),/get_lun
   'l': openr,lunr_rl,string(disk+'runlogs/lab/'+occul),/get_lun
endcase

readf,lunr_rl,nhl,ncol   ; title line
ihl=1
while(ihl lt nhl) do begin
readf,lunr_rl,rlrecord   ; skip title line
ihl=ihl+1
endwhile

ll_year=strpos(rlrecord,'Year')
nspec=0l
print, 'locating SPT files.....'
;  Loop over spectra in runlog
while(eof(lunr_rl) eq 0) do begin
   readf,lunr_rl,rlrecord
   prefix=strmid(rlrecord,0,1)
   if ( prefix ne ":" ) then begin
      if ( prefix eq " " or prefix eq "+" or prefix eq "-" ) then rlrecord=strmid(rlrecord,1,ll_year-2)
      spectrum=strtrim(strmid(rlrecord,0,ll_year-2),2)
      print,spectrum
      sptpath(nspec)=disk+'spt/z'+spectrum
      nspec=nspec+1
   endif
;   print,nspec-1,sptpath(nspec-1)
endwhile
print,nspec,' SPT files located.'
if nspec gt mspec then print, 'Increase parameter MSPEC=',mspec
close,lunr_rl
free_lun,lunr_rl

loadct,39
kspec=-1
kskip=1
two:
  kspec=kspec+kskip
  if kspec lt 0 then begin
      print,format='(" Off beginning of list. Gone back to first spectrum.")'
      kspec=0
      kskip=+1
  endif
  if kspec gt nspec-1 then begin
      print,format='(" Off end of list. Gone back to last spectrum.")'
      kspec=nspec-1
      kskip=-1
  endif

  xunit=0 ; cm-1 (default)
; Open next SPT file
  sptpathlen=strlen(findfile(sptpath(kspec)))
  if sptpathlen(0) le 0  then goto,two
  openr,lunr_spt,sptpath(kspec), /get_lun
  on_ioerror,closeunit
  readf,lunr_spt,nhl,ncol
  ntgas=ncol-4  ; the first 4 columns are Freq, Tm, Tc, Cont
;  print,nhl,ncol,sptpath(kspec),ntgas
  if(ntgas gt 2) then begin
     readf,lunr_spt,format='(2f14.6,i7,f8.4,3f8.3,1x,f9.6,e10.3,f7.4,f8.4,3(f7.3,1x,i3))',$
     fmin,fmax,npoints,effres,asza,zobs,tang,rmsocl,peff,frac,xzo,colmant,colexp,vsfmant,vsfexp,errmant,errexp
;     print,fmin,fmax,npoints,effres,asza,zobs,tang,rmsocl,peff,frac,colmant,colexp,errmant,errexp
;     expont=max([colexp,errexp])
;     burden=colmant*vsfmant*(10^(colexp-expont))
;     berr=colmant*vsfmant*errmant*(10^(errexp-expont))
     burden=colmant*vsfmant*10^(vsfexp)
     berr=colmant*errmant*10^(errexp)
;     print,'column=',burden,berr,colexp
  endif else begin
     readf,lunr_spt,format='(2f14.6,i7,f8.4,3f8.3,1x,f9.6,e10.3,f7.4,f8.4)',$
     fmin,fmax,npoints,effres,asza,zobs,tang,rmsocl,peff,frac,xzo
     burden=0.0
     berr=0.0
     colexp=0
  endelse
;  print,'sptpath,npoints=',sptpath(kspec),npoints
  if npoints le 0 then goto,two
  readf,lunr_spt,headers
;  print,'headers=',headers
;  headers=headers+'other'
  if kskip ne 0 then begin
  if same_window eq 1 then begin
    vbar=0.5*(fmin+fmax)
    if abs(vbar-vbar_was) gt 0.5*(fmax-fmin) then begin
      close,lunr_spt
      free_lun,lunr_spt
      goto,two
    endif
  endif
  endif
  vbar_was=0.5*(fmin+fmax)
closeunit:
  print,'After closeunit'
  print,"defining datarray: ",4+ntgas,npoints
  datarray=dblarr(4+ntgas,npoints)
;  print,"defined datarray: ",npoints
  readf,lunr_spt,datarray

  datarray(1,*)=(datarray(1,*)/datarray(3,*)-xzo)/(1.0-xzo) ; (obsrvd(k)/cont(k)-xzo)/(1-xzo),
  datarray(2,*)=(datarray(2,*)/datarray(3,*)-xzo)/(1.0-xzo) ; (calcul(k)/cont(k)-xzo)/(1-xzo)
  close,lunr_spt
  free_lun,lunr_spt
  print," closed lunr_spt: ",npoints
  if npoints le 0 then goto,two
;
;; The following loop multiplies the various target transmittances
;; by the "other" transmittances. This is useful to bring all the
;; contributions to the same tramsmittance level when "other" is
;; a strong continuum absorption.
;  for icol=3,2+ntgas do begin
;    datarray(icol,*) = datarray(icol,*)*datarray(2+ntgas,*)
;  endfor
;
  if abs(burden) gt 99.999  then begin
    burden=burden/10
    berr=berr/10
    colexp=colexp+1
  endif

  indszp2=sptpath(kspec).Indexof('/z')+2
  sptpathlen=strlen(sptpath(kspec))-indszp2
  text(ntgas)='  '+string(strmid(sptpath(kspec),indszp2,min([37,sptpathlen])))+$
  string(format='("   !7w!6 =",(f6.2),"!9%")',asza)+$
  string(format='("   !6Z!dT!n =",(f7.2),"km")',tang)+$
  string(format='("   !7r!6!drms!n =",(f7.4),"%")',100*rmsocl)
;  string(format='("  !9i!6dz=",(f6.3))',burden)+$
;  '!9+!6'+string(format='(f6.3)',berr)+$
;  string(format='("x10!u",(i2))',colexp)

;  find the highest and lowest data point for each category
  yloarr=fltarr(ntgas+4)
  yhiarr=fltarr(ntgas+4)
;  print,'npoints,ntgas=',npoints,ntgas
  for tt=0,ntgas+3 do begin
    yloarr(tt)=min(datarray(tt,0:npoints-1))
    yhiarr(tt)=max(datarray(tt,0:npoints-1))
;    print,'tt,ylo,yhi= ',tt,yloarr(tt),yhiarr(tt)
  endfor

while idn(0) lt 24 do begin

  if( ktgas eq ntgas+1) then begin ; Continuum plotted
     ylo=min(yloarr(3))
     yhi=max(yhiarr(3))
     print,'0,ktgas,ntgas,ylo,yhi: ',ktgas,ntgas,ylo,yhi
  endif else if( ktgas eq ntgas) then begin ; All target gases plotted
     print,'1,ktgas,ntgas,ylo,yhi: ',ktgas,ntgas,ylo,yhi
     ylo=min([yloarr(1),yloarr(2),yloarr(4:ntgas+3)])
     yhi=max([yhiarr(1),yhiarr(2),yhiarr(4:ntgas+3)])
  endif else if (ktgas ge 0) then begin ; One target gases plotted
     if ( ktgas gt ntgas ) then  ktgas=0
     print,'2,ktgas,ntgas,ylo,yhi: ',ktgas,ntgas,ylo,yhi
     ylo=min([yloarr(1),yloarr(2),yloarr(ktgas+4)])
     yhi=max([yhiarr(1),yhiarr(2),yhiarr(ktgas+4)])
  endif else begin         ;  Zero target gases plotted (just Tm & Tc)
;     if ( ktgas gt ntgas ) then  ktgas=0
     ylo=min([yloarr(1),yloarr(2)])
     yhi=max([yhiarr(1),yhiarr(2)])
     print,'3,ktgas,ntgas,ylo,yhi: ',ktgas,ntgas,ylo,yhi
  endelse
  if rlim eq "a" then  rmax=max(abs(datarray(2,0:npoints-1)-datarray(1,0:npoints-1)))
  if xlolim eq "a" then  xmin=fmin
  if xuplim eq "a" then  xmax=fmax
  if ylolim eq "a" then  ymin=ylo
  if yuplim eq "a" then  ymax=yhi
  !x.range=[1.001*xmin-0.001*xmax,1.001*xmax-0.001*xmin]
  !y.range=[1.001*ymin-0.001*ymax,1.001*ymax-0.001*ymin]
;
;  Write gas names to right of main panel.
  kpoints=npoints*(xmax-fmin)/(fmax-fmin)
gases = strsplit(headers,count=nss,/extract)
;print, ntgas,gases
;  substr, headers, gases, mss, nss
  for i=0, ntgas-1 do begin
     if mtxt(i) eq 0 then begin
        xtxt(i)=xmax+0.004*(xmax-xmin)
        ytxt(i)=ymax-(i+0.5)*(ymax-ymin)/15
        if mtxt(i) eq 0 then text(i)=gases(i+4)  ; first 3 columns are Freq, Tm, Tc, Cont
     endif
  endfor
;
  time=strmid(systime(0),8,2)+'-'+strmid(systime(0),4,3)+'-'+$
  strmid(systime(0),20,4)+'/'+strmid(systime(0),11,8)
;  if strlen(text(ntgas+1)) le 0 then text(ntgas+1)=filenm+'  '+time
 print,'ylolim,ymin=',ylolim,ymin
 print,'yuplim,ymax=',yuplim,ymax
  pplot,xmin,xmax,ymin,ymax,datarray,ntgas,npoints,nbox,rmax,ktgas,xzo,xunit,xtxt,ytxt,text,ntgas+2+ntxt
;
  event=widget_event(b1)
  idn=where(butn eq event.id)
  case idn(0) of
    0: begin     ; N-2
          kskip=-2
          goto,two
       end
    1: begin     ; N-1
          kskip=-1
          goto,two
       end
    2: begin     ; replot
          kskip=+0
          goto,two
       end
    3: begin     ; N+1
          kskip=+1
          goto,two
       end
    4: begin     ; N+2
          kskip=+2
          goto,two
       end
    5: begin     ; skip
          print,format='("currently at file ",(i4)," of ",(i4))',kspec,nspec
          print,format='($,"how many files to skip?")
          read,kskip
          goto,two
       end
    6: begin     ; residual
          nbox=3-nbox  ;  toggles between 1 & 2
       end
    7: begin     ; R_max
          print,format='($,"Max R-value")'
          read,rmax
          rlim='m'
       end
    8: begin     ; R_auto 
          rlim='a'
       end
    9: begin     ; X_min
          print,format='($,"Min X-value")'
          read,xmin
          xlolim='m'
       end
   10: begin     ; X_max
          print,format='($,"Max X-value")'
          read,xmax
          xuplim='m'
       end
   11: begin     ; X_auto 
          xuplim='a'
          xlolim='a'
       end
   12: begin     ; X_left
          dx=xmax-xmin
          xmax=xmin
          xmin=xmin-dx
          xlolim='m'
          xuplim='m'
       end
   13: begin     ; X_right
          dx=xmax-xmin
          xmin=xmax
          xmax=xmax+dx
          xlolim='m'
          xuplim='m'
       end
   14: begin     ; X_zoom
          temp=0.75*xmax+0.25*xmin
          xmin=0.25*xmax+0.75*xmin
          xmax=temp
          xlolim='m'
          xuplim='m'
       end
   15: begin     ; X_anti-zoom
          temp=1.5*xmax-0.5*xmin
          xmin=1.5*xmin-0.5*xmax
          xmax=temp
          xlolim='m'
          xuplim='m'
       end
   16: begin     ; Y_min
          print,format='($,"Min Y-value")'
          read,ymin
          ylolim='m'
       end
   17: begin     ; Y_max
          print,format='($,"Max Y-value")'
          read,ymax
          yuplim='m'
       end
   18: begin     ; Y_auto 
          yuplim='a'
          ylolim='a'
       end
   19: begin     ; target gases
;          ktgas = ( (ktgas+2) mod (ntgas+2) ) -1
          ktgas = ( (ktgas+2) mod (ntgas+3) ) -1
       end
   20: begin     ; Text
          for itxt=0,ntgas+2+ntxt-1 do print,format='(i2,2e15.6,1x,a)',itxt,xtxt(itxt),ytxt(itxt),text(itxt)
         print,' Enter #N, X, Y, Text'
         read, itxt,xt,yt,sss
         text(itxt)=strmid(sss,1,999)  ; skip  SSS(0) which is the delimiter
         xtxt(itxt)=xt
         ytxt(itxt)=yt
         ntxt=max([itxt+1-ntgas-2,ntxt])
         mtxt(itxt)=1
       end
   21: begin     ; same window
          same_window=1-same_window   ; toggles between 0 and 1
          vbar_was=0.5*(fmin+fmax)
       end
   22: begin     ; cm-1/nm
          print,fmin,fmax,xmin,xmax,xlolim,xuplim
          datarray(0,*)=10000000/datarray(0,*)
          xunit=1-xunit

          temp=10000000/fmax
          fmax=10000000/fmin
          fmin=temp

          temp=10000000/xmax
          xmax=10000000/xmin
          xmin=temp

          atemp=xlolim
          xlolim=xuplim
          xuplim=atemp
          print,fmin,fmax,xmin,xmax,xlolim,xuplim
       end
   23: begin     ; hardcopy
          set_plot,'ps'
          filenm='pfit_'+string(format='(i3.3)',psn)+$
           '.'+strlowcase(!d.name)
;          device,/landscape,/color,font_size=7,filename=filenm
          device,/portrait,/color,font_size=7,filename=filenm
          !x.range=[1.001*xmin-0.001*xmax,1.001*xmax-0.001*xmin]
          !y.range=[1.001*ymin-0.001*ymax,1.001*ymax-0.001*ymin]
          pplot,xmin,xmax,ymin,ymax,datarray,ntgas,npoints,nbox,rmax,ktgas,xzo,xunit,xtxt,ytxt,text,ntgas+2+ntxt
          device,/close
          set_plot,idl_device
          psn=psn+1
       end
   24: print, ' Quitting plot program'
  endcase
endwhile

widget_control,b1,/destroy
wdelete,0
end
