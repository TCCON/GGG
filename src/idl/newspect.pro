pro newspect

color=string('f')
idl_device=string('x')
occul=string(14)
resp=string('n')
runlogstr=strarr(80000)
spectrum=string('                                   ')
termtp=string('xwin')
names=''

apo=0
ipart=0
npart=0
kspec=0l
ifirst=0l
ilast=0l
graw=0d
path=string(' ')
dwas=' '
hdcp=99
intrp=1
lambdaflag=0
oplotflag=0
ll=0
llim=string('f')
nocl=1
ttyp=0
possp=0l
ulim=string('f')
ifset=0
version=string(' newspect    Version 2.10     2018-04-17    gct ')

!p.charsize=1.5
!p.thick=3
!x.thick=2
!y.thick=2
!p.charthick=2

!x.minor=10
!y.style=1
!y.minor=5
!ytitle='!17Signal (arbitrary units)'

print,version
;termin,ttyp,termtp,color,idl_device
set_plot,idl_device
window,0,title='newspect.pro', xsize=1100,ysize=512,colors=2
  xmenu,base=base,['prev (-2)','prev (-1)','next (+1)','next (+2)','skip',$
  'yscale upper','yscale lower','yscale auto','frequency','zoom','antizoom',$
  'cursor','frame right','frame left','freq+4','atlas frame','apodize',$
  'interpolate','w.length','oplot','hardcopy','quit'],buttons=buttn
widget_control,/realize,base

rtn1:
on_ioerror,rtn1
print,format='($,"enter runlog")'
read,occul
if occul eq 'q' then goto,rtn5
;idot=strpos(occul,'.')
;case strmid(occul,idot+1,1) of
ix=strlen(occul)
case strmid(occul,ix-3,1) of
   'g': openr,unit,string('$GGGPATH/runlogs/gnd/',occul),/get_lun
   'b': openr,unit,string('$GGGPATH/runlogs/bal/',occul),/get_lun
   'o': openr,unit,string('$GGGPATH/runlogs/orb/',occul),/get_lun
   'l': openr,unit,string('$GGGPATH/runlogs/lab/',occul),/get_lun
endcase

;print,format='($,"central frequency (cm-1)")'
;read,frqcen
;print,format='($,"window width (cm-1)")'
;read,wid

on_ioerror,rtn5

readf,unit,nhl,ncol   ; title line
ihl=1
while(ihl lt nhl) do begin
readf,unit,names   ; skip title line
  print,ihl,names
ihl=ihl+1
endwhile

;readf,unit,names
;readf,unit,names
;readf,unit,names
;print,strpos(names,'Year')
ll_year=strpos(names,'Year')
rlformat=string('(1x,a',ll_year-2,',34x,2f8.3,f7.4,f8.3,f7.3,f7.2,3f6.4,2i9,f15.11,i9,i3,39x,f8.2)')
;print,rlformat
j=0l
while(eof(unit) eq 0) do begin
  readf,unit,names
  if strmid(names,0,7) eq 'Format=' then rlformat=strmid(names,7)
  runlogstr(j)=names
  j=j+1
endwhile
nspe=j
close,unit
free_lun,unit
;print,format='($,"spectrum number (-1 to quit)")'
;read,kspec
print,'nspe=',nspe
if kspec lt 0 or kspec gt nspe-1 then goto,rtn5
kspec=kspec-1
kskip=1l

rtn1b:
kskip=kskip/abs(kskip)   ; +1 or -1

rtn2:
kspec=kspec+kskip
if kspec lt 0 then begin
   print,format='(" Off beginning of list. ")'
   kspec=0
   kskip=1
endif
if kspec gt nspe-1 then begin
   print,format='(" Off end of list. ")'
   kspec=nspe-1
   kskip=-1
endif
print,'kspec,kskip=',kspec,kskip
reads,runlogstr(kspec),format=rlformat,$
spectrum,zobs,asza,zenoff,azim,osds,opd,fovi,fovo,amal,ifirst,ilast,graw,possp,bpw,pout

if (ifset le 0) then begin  ; No freq range defined, so plot all.
  frqcen=graw*(ifirst+ilast)/2
  wid=graw*(ilast-ifirst)
endif

spectrum=strtrim(spectrum)
;print, spectrum,zobs, asza+zenoff, (3389.5+zobs)*sin(3.14159265*(asza+zenoff)/180)-3389.5
;graw=graw*(double(1.)+(fovi^2+amal^2)*6.25e-02)
if ((round(frqcen-0.5*wid) gt ilast*graw) or $
  (long(frqcen+0.5*wid) lt ifirst*graw)) then begin
  print,format='(a22," only encompasses ",f9.1," to ",f9.1, " cm-1")',$
    spectrum,ifirst*graw,ilast*graw
endif
ffile,spectrum,partition,ipart,npart,path
;ffile,spectrum,currpart,path,flag
if strlen(path) eq 0 then begin
  print,format='(a22," not found in data partitions")',spectrum
  print,format='("reading next spectrum")'
  goto,rtn2
endif

rtn3:
readsp,path,frqcen-(0.5*wid),frqcen+(0.5*wid),opd,ifirst,ilast,$
  graw,possp,bpw,apo,intrp,freq,yobs,nobs,iflag
if iflag ne 0 then begin
  print,format='(" Problem reading ",a22)',spectrum
  print,format='("trying next spectrum")'
  goto,rtn2
endif

rtn4a:
if llim eq 'f' then mintrans=min(yobs)
if ulim eq 'f' then maxtrans=max(yobs)

rtn4b:
!y.range=[mintrans,maxtrans]
!mtitle=string(spectrum)+$
  string(format='("  !7v!6:",(f6.2),"!9% ")',(asza+zenoff))+$
  string(format='("  !6pa:",(f7.2),"hPa")',pout)+$
  string(format='("  apo:",(i1))',apo)+$
  string(format='("  intrp:",(i2))',intrp)
!x.style=1
!x.ticklen=-.03
!y.ticklen=-.03
if lambdaflag eq 0 then begin
  !x.range=[frqcen-0.5*wid,frqcen+0.5*wid]
  !xtitle='!17Frequency (cm!u-1!n)'
;  plot,freq,yobs,linestyle=0,noclip=nocl
if oplotflag eq 0 then  plot,freq,yobs,linestyle=0,noclip=nocl else oplot,freq,yobs,linestyle=0,noclip=nocl
endif else begin
  !x.range=[1.e+04/(frqcen-0.5*wid),1.e+04/(frqcen+0.5*wid)]
  !xtitle='Wave length (!7l!6m)'
;  plot,(1.e+04/freq),yobs,linestyle=0,noclip=nocl
if oplotflag eq 0 then  plot,(1.e+04/freq),yobs,linestyle=0,noclip=nocl else oplot,(1.e+04/freq),yobs,linestyle=0,noclip=nocl
endelse

print,oplotflag,lambdaflag

;for j=0,nobs-1 do begin
;    print,j,freq(j),yobs(j)
;endfor 

time=strmid(systime(0),8,2)+'-'+strmid(systime(0),4,3)+'-'+$
  strmid(systime(0),20,4)+'/'+strmid(systime(0),11,8)
  xyouts,930,10,time,charsize=0.8,/device
event=widget_event(base)
ll=where(buttn eq event.id)
case ll(0) of
;   0: begin          ; runlog
;        dwas=strmid(occul,strlen(occul)-2,2)
;        goto,rtn1
;      end
   0: begin          ; previous -2
        kskip=-2
        goto,rtn2
      end
   1: begin          ; previous
        kskip=-1
        goto,rtn2
      end
   2: begin          ; next
        kskip=+1
        goto,rtn2
      end
   3: begin          ; next +2
        kskip=+2
        goto,rtn2
      end
   4: begin          ; skip
        print,format='("currently at spectrum #:",(i5)," /",(i5))',kspec,nspe
        print,format='($,"how many spectra to skip?")'
        read,kskip
        goto,rtn2
      end
   5: begin          ; yscale upper
        print,format='($,"transmission maximum for window")'
        read,maxtrans
        ulim=string('t')
        goto,rtn4b
      end
   6: begin          ; yscale lower
        print,format='($,"transmission minimum for window")'
        read,mintrans
        llim=string('t')
        goto,rtn4b
      end
   7: begin          ; yscale auto
        llim=string('f')
        ulim=string('f')
        goto,rtn4a
      end
   8: begin          ; frequency
        print,format='($,"central frequency (cm-1)")'
        read,frqcen
        print,format='($,"window width (cm-1)")'
        read,wid
        ifset=1
        goto,rtn3
      end
   9: begin          ; zoom
        wid=wid*0.5
        ifset=1
        goto,rtn3
      end
   10: begin          ; antizoom
        wid=wid*2.0
        ifset=1
        goto,rtn3
      end
   11: begin          ; cursor
        cursor,x,y,/data,/wait
        frqcen=x
        print,format='("central frequency (cm-1) ",(f9.4))',frqcen
        print,format='($,"window width (cm-1)")'
        read,wid
        goto,rtn3
      end
  12: begin          ; frame right
        frqcen=frqcen+wid
        ifset=1
        goto,rtn3
      end
  13: begin          ; frame left
        frqcen=frqcen-wid
        ifset=1
        goto,rtn3
      end
  14: begin          ; freq+4
        frqcen=frqcen+4.0
        ifset=1
        goto,rtn3
      end
  15: begin          ; atlas frame
        slf,yobs,nobs,y0,grad
        for i=0,nobs-1  do begin
           yobs(i)=yobs(i)/(y0+grad*i)
        endfor
        yobs=[yobs,rotate(10*yobs-9,2)]
        freq=[freq,rotate(freq,2)]
        goto,rtn4a
      end
  16: begin          ; apodise
        print,format='("apodisation currently set at: ",i3)',apo
        print,format='($,"apodize [0-3]")'
        read,apo
        goto,rtn3
      end
  17: begin          ; interpolate
        print,format='("interpolation currently set at: ",i3)',intrp
        intrp=0
        while intrp le 0 do begin
           print,format='($,"interpolate [1-10]")'
           read,intrp
        endwhile
        goto,rtn3
      end
  18: begin          ; w.length
        case lambdaflag of
          0: lambdaflag=1
          1: lambdaflag=0
        endcase
        goto,rtn4b
      end
  19: begin          ; oplotflag
        case oplotflag of
          0: oplotflag=1
          1: oplotflag=0
        endcase
        goto,rtn4b
      end
  20: begin          ; hardcopy
        set_plot,'ps'
        spawn,'ls -l *.ps',pslist
        num_add=max(where(pslist))
        if num_add lt 0 then psn=hdcp+1
        if num_add ge 0 then psn=hdcp+num_add+2
        filenm='newspect_'+string(format='((i3))',psn)+$
          '.'+strlowcase(!d.name)
        device,/portrait,filename=filenm
        print,format='("programme writing postscript file")'
        print,format='("please wait for screen update")'
        if lambdaflag eq 0 then begin
          !x.range=[frqcen-0.5*wid,frqcen+0.5*wid]
          !xtitle='!17Frequency (cm!u-1!n)'
;          plot,freq,yobs,linestyle=0,noclip=nocl
          plot,freq,yobs,linestyle=0
        endif else begin
          !x.range=[1.e+04/(frqcen-0.5*wid),1.e+04/(frqcen+0.5*wid)]
          !xtitle='Wave length (!7l!6m)'
;          plot,(1.e+04/freq),yobs,linestyle=0,noclip=nocl
          plot,(1.e+04/freq),yobs,linestyle=0
        endelse
        time=strmid(systime(0),8,2)+'-'+strmid(systime(0),4,3)+'-'+$
          strmid(systime(0),20,4)+'/'+strmid(systime(0),11,8)
        xyouts,930,10,time,charsize=0.8,/device
        device,/close
;        spawn,'lp '+filenm
        set_plot,idl_device
        goto,rtn4a
      end
  21: begin          ; quit
        goto,rtn5
      end
endcase

rtn5:
print,'end of idl programme'
wdelete,0
widget_control,base,/destroy
end
