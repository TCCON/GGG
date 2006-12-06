      pro lnbc,label,nc
; Returns position of Last Non-Blank Character in string LABEL
; Recognized "blank" characters include:
;       Null              ASCII # 00
;       Horizontal Tab    ASCII # 09
;       Space             ASCII # 32
;       Comma             ASCII # 44
; NC=-1 indicates that the entire string was blank^
      for nc=strlen(label)-1,0,-1  do begin
         cc=strmid(label,nc,1)
         if (    cc ne string(00B)
            and  cc ne string(09B)
            and  cc ne string(32B)
            and  cc ne string(44B) ) then return
      endfor
      nc=-1
      return       ;  Abnormal return: No characters found
      end
