pro newigram

fname=string(replicate(32b,12))
idl_device='x'
cmd=' '
color='f'
chg=0
cin=0
cx=0
gmt=' '
hdcp=99
length=' '
location=' '
pdat=create_struct('ipart',0,'npart',0)
prevcmd=' '
resp='y'
termtp='xwin'
ttyp=0
type=' '
version='newig  v.1.0.0, 05/05/98'

print,version
print,format='($,"interferogram")
read,fname
newffile,fname,pdat,glenfiddich,flag
openr,unit,glenfiddich,/get_lun
igram=fstat(unit)
xsize=long(igram.size/256)
print,format='("number of frames ",(i5))',xsize
case (igram.size mod 128) of
     0: xvalue=intarr(xsize*128)
  else: xvalue=intarr((xsize+1)*128)
endcase
readu,unit,xvalue
;xvalue=swap_endian(xvalue)
yvalue=intarr(128)
termin,ttyp,termtp,color,idl_device
set_plot,idl_device
repeat begin
  print,format='($,"enter au,bo,cx,hc,hd,hg,in,qu,rd,frm,yhg,yin")
  read,cmd
  case cmd of
     'au': begin
             if cin eq 1 then begin
               print,format='($,"set insb yscale automatically [y/n]")'
               read,resp
               if resp eq 'y' then cin=0
             endif
             if chg eq 1 then begin
               print,format='($,"set hgcd yscale automatically [y/n]")'
               read,resp
               if resp eq 'y' then chg=0
             endif
           end
     'bo': begin
             !p.multi=[0,0,2,0,0]
             !p.title=$
               string(fname,'    ',gmt,'    ',length,' ',type,' ',location)
             !x.range=[frame(0),frame(1)]
             if cx eq 1 then !x.range=[xmin,xmax]
             !x.style=1
             !x.tickformat='(i8)'
             !x.title='HgCdTe Channel'
             !y.range=[hgmin,hgmax]
             !y.style=16
             if chg eq 1 then !y.range=[yhgmin,yhgmax]
             plot,xarray(0:k(0)-1)/60.,hgcd(0:k(0)-1),linestyle=0
             !p.title=''
             !x.title='InSb Channel'
             !y.range=[inmin,inmax]
             if cin eq 1 then !y.range=[yinmin,yinmax]
             plot,xarray(0:k(1)-1)/60.,insb(0:k(1)-1),linestyle=0
             !p.multi=0
           end
     'cx': begin
             print,format='($,"enter xmin and xmax")'
             read,xmin,xmax
             cx=1
           end
     'hc': begin
             set_plot,'ps'
             spawn,'ls -l idl_*.ps',pslist
             num_add=max(where(pslist))
             if num_add lt 0 then psn=hdcp+1
             if num_add ge 0 then psn=hdcp+num_add+2
             filenm='idl_'+string(format='(i3)',psn)+$
               '.'+strlowcase(!d.name)
             device,/landscape,font_size=11,filename=filenm
             !p.title=$
               string(fname,'    ',gmt,'    ',length,' ',type,' ',location)
             case prevcmd of
                'bo': begin
                        !p.multi=[0,0,2,0,0]
                        !x.range=[frame(0),frame(1)]
                        if cx eq 1 then !x.range=[xmin,xmax]
                        !x.style=1
                        !x.tickformat='(i8)'
                        !x.title='HgCdTe Frames'
                        !y.range=[hgmin,hgmax]
                        if chg eq 1 then !y.range=[yhgmin,yhgmax]
                        !y.style=16
                        plot,xarray(0:k(0)-1)/60.,hgcd(0:k(0)-1),psym=3
                        !p.title=''
                        !x.title='InSb Frames'
                        !y.range=[inmin,inmax]
                        if cin eq 1 then !y.range=[yinmin,yinmax]
                        plot,xarray(0:k(1)-1)/60.,insb(0:k(1)-1),psym=3
                        !p.multi=0
                      end
                'hg': begin
                        !x.range=[frame(0),frame(1)]
                        if cx eq 1 then !x.range=[xmin,xmax]
                        !x.style=1
                        !x.tickformat='(i8)'
                        !x.title='HgCdTe Frames'
                        !y.range=[hgmin,hgmax]
                        if chg eq 1 then !y.range=[yhgmin,yhgmax]
                        !y.style=16
                        plot,xarray(0:k(0)-1)/60.,hgcd(0:k(0)-1),linestyle=0
                      end
                'in': begin
                        !x.range=[frame(0),frame(1)]
                        if cx eq 1 then !x.range=[xmin,xmax]
                        !x.tickformat='(i8)'
                        !x.title='InSb Frames'
                        if cin ne 1 then !y.range=[inmin,inmax]
                        !y.style=16
                        plot,xarray(0:k(1)-1)/60.,insb(0:k(1)-1),linestyle=0
                      end
                else: print,'prevcmd not a plot command...'
             endcase
             device,/close
             spawn,'lp '+filenm
             set_plot,idl_device
           end
     'hd': begin
             for z=1,10 do begin
               yvalue=xvalue(128*z:127+128*z)
               ;print,format='(129z5.4)',z,yvalue
               case (yvalue(121) and '3c00'x) of
                 '2000'x: begin
                            header_frame=header_frame+1
                            gmt=string((ishft(yvalue(0),-4) and '000f'x))+$
                                  string((yvalue(0) and '000f'x))+$
                                  string((ishft(yvalue(1),-4) and '000f'x))+$
                                  string((yvalue(1) and '000f'x))+$
                                  string((ishft(yvalue(2),-4) and '000f'x))+$
                                  string((yvalue(2) and '000f'x))+$
                                  string((ishft(yvalue(3),-4) and '000f'x))+$
                                  string((yvalue(3) and '000f'x))+$
                                  string((ishft(yvalue(4),-4) and '000f'x))+$
                                  string((yvalue(4) and '000f'x))+$
                                  string((ishft(yvalue(5),-4) and '000f'x))+$
                                  string((yvalue(5) and '000f'x))
                            gmt=strcompress(gmt,/remove_all)
                            a=(yvalue(6) and '0010'x)
                            case a of
                                 0: type='normal'
                              else: type='period'
                            endcase
                            a=(yvalue(8) and '0002'x)
                            case a of
                                 0: type=type+' reverse run'
                              else: type=type+' forward run'
                            endcase
                            a=(yvalue(6) and '0007'x)
                            case a of
                                 0: length='  4 cm.'
                                 1: length=' 30 cm.'
                                 2: length=' 60 cm.'
                                 3: length='120 cm.'
                                 4: length='180 cm.'
                              else: length='  0 cm.'
                            endcase
                            a=ishft((yvalue(6) and '00e0'x),-5)
                            case a of
                                 0: location='@ AIR'
                                 1: location='@ JPL'
                                 2: location='@ MES'
                                 3: location='@ TMO'
                                 4: location='@ AMS'
                                 5: location='@ PAL'
                                 6: location='@ CHI'
                                 7: location='@ ANT'
                              else: location=' UNKN'
                            endcase
                            goto,tmp
                          end
                 else:
               endcase
             endfor
           end
     'hg': begin
             !p.title=$
               string(fname,'    ',gmt,'    ',length,' ',type,' ',location)
             !x.range=[frame(0),frame(1)]
             if cx eq 1 then !x.range=[xmin,xmax]
             !x.style=1
             !x.tickformat='(i8)'
             !x.title='HgCdTe Frames'
             !y.range=[hgmin,hgmax]
             if chg eq 1 then !y.range=[yhgmin,yhgmax]
             !y.style=16
             plot,xarray(0:k(0)-1)/60.,hgcd(0:k(0)-1),psym=3
           end
     'in': begin
             !p.title=$
               string(fname,'    ',gmt,'    ',length,' ',type,' ',location)
             !x.range=[frame(0),frame(1)]
             if cx eq 1 then !x.range=[xmin,xmax]
             !x.style=1
             !x.tickformat='(i8)'
             !x.title='InSb Frames'
             !y.range=[inmin,inmax]
             if cin eq 1 then !y.range=[yinmin,yinmax]
             !y.style=16
             plot,xarray(0:k(1)-1)/60.,insb(0:k(1)-1),psym=3
           end
     'qu': begin
             print,'later. bye...'
             tmp1:
             close,unit
             free_lun,unit
           end
     'rd': begin
             on_ioerror,tmp1
             for z=(frame(0)-1),(frame(1)-1) do begin
               yvalue=xvalue(128*z:127+128*z)
               ;print,format='(129z5.4)',z,yvalue
               q=z-(frame(0)-1)
               ;print,$
               ;  format='("i(121), i(126), i(127), i(121) & x03c00 ",4z)',$
               ;  yvalue(121),yvalue(126),yvalue(127),(yvalue(121) and '3c00'x)
               case (yvalue(121) and '3c00'x) of
                 '0400'x: begin
                            here:
                            filler_frame=filler_frame+1
                            insb(60*q:60*q+59)='0800'x
                            hgcd(60*q:60*q+59)='0800'x
                            k(0)=60*q+59
                            k(1)=60*q+59
                          end  ; end of case '0400'x
                 '0800'x: begin
                            science_frame=science_frame+1
                            for i=0, 119 do begin
                              case (i and '0001'x) of
                                '0'x: begin
                                        hgcd(k(0))=yvalue(i)
                                        ; value=(yvalue(i) and '0fff'x)-'0800'x
                                        ; range=ishft((yvalue(i) and '7000'x),$
                                        ;   -12)
                                        ; hgcd(k(0))=ishft(long(value),range)
                                        k(0)=k(0)+1
                                      end
                                '1'x: begin
                                        insb(k(1))=yvalue(i)
                                        ;value=(yvalue(i) and '0fff'x)-'0800'x
                                        ;range=ishft((yvalue(i) and '7000'x),$
                                        ;  -12)
                                        ;insb(k(1))=ishft(long(value),range)
                                        k(1)=k(1)+1
                                      end
                              endcase
                            endfor
                          end  ; end of case '0800'x
                 '1000'x: begin
                            ; print,format='("1000x option: trailer frame")'
                            trailer_frame=trailer_frame+1
                            insb(60*q:60*q+59)='0800'x
                            hgcd(60*q:60*q+59)='0800'x
                            k(0)=60*q+59
                            k(1)=60*q+59
                          end  ; end of case '1000'x
                 '2000'x: begin
                            header_frame=header_frame+1
                            insb(60*q:60*q+59)='0800'x
                            hgcd(60*q:60*q+59)='0800'x
                            k(0)=60*q+59
                            k(1)=60*q+59
                          end  ; end of case '2000'x
                    else: begin
                            goto, here
                          end
               endcase
             endfor
             mantissa=n_elements(hgcd)
             range=n_elements(hgcd)
             mantissa=(hgcd and '0fff'x)-'0800'x
             range=ishft((hgcd and '7000'x),-12)
             hgcd=ishft(long(mantissa),range)
             mantissa=(insb and '0fff'x)-'0800'x
             range=ishft((insb and '7000'x),-12)
             insb=ishft(long(mantissa),range)
             print,format='("number of header frames  ",(i6))',header_frame
             print,format='("number of trailer frames ",(i6))',trailer_frame
             print,format='("number of filler frames  ",(i6))',filler_frame
             print,format='("number of science frames ",(i6))',science_frame
             inmin=min(insb)
             inmax=max(insb)
             hgmin=min(hgcd)
             hgmax=max(hgcd)
             ;print,inmin,inmax,hgmin,hgmax
           end
    'frm': begin
             k=[long(0),long(0)]
             header_frame=0
             trailer_frame=0
             filler_frame=0
             science_frame=0
             print,format='($,"enter start and stop frame")'
             read,m1,m2
             frame=[long(m1),long(m2)]
             delta=(frame(1)-frame(0))/10.
             insb=lonarr(delta*600+60)
             hgcd=lonarr(delta*600+60)
             xarray=fltarr(long(delta*600+60))
             for i=long(0),delta*600+59 do xarray(i)=float(i+60*(frame(0)-1))
           end
    'yin': begin
             print,format='($,"enter insb ymin and ymax")'
             read,yinmin,yinmax
             cin=1
           end
    'yhg': begin
             print,format='($,"enter hgcdte ymin and ymax")'
             read,yhgmin,yhgmax
             chg=1
           end
     else: print,'command not recognized...'
  endcase
  tmp:
  prevcmd=cmd
endrep until cmd eq 'qu'
end
