pro getendian,iend

d=131071l
iend=intarr(1)

openw,unit,'t',/get_lun
writeu,unit,d
point_lun,unit,0
readu,unit,iend
close,unit
free_lun,unit

spawn,'/bin/rm -f t'
end
