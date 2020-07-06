pro getendian,iend
;  Determines the "endian-ness" of the host computer
;     iend=+1 on big-endian CPU (e.g. Sun)
;     iend=-1 on little-endian CPU (e.g. PC)
; 131071 = 2^17-1, is the number whose first 2 bytes = +1 and whose second 2 bytes = -1.

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
