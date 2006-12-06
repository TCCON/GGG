      pro mmm,ipanel,ncol,nrow,dir,frx,fry,x0,y0,x1,y1
;
; INPUTS:
;        IPANEL  Panel Number (from 1 to NCOL*NROW)
;        NCOL  Numer of columns
;        NROW  Number of ROWs
;        DIR   Flag (0 or 1) indicating the direction of panel increments
;        FRX   Fraction of x-space to be used for data window
;        FRY   Fraction of y-space to be used for data window
;
; OUTPUTS
;        X0, Y0  Co-ordinate of origin of data window
;        X1, Y1  Co-ordinate of upper right of data window
;
; Note that a margin of width 0.1 is alwas left along the
; left and bottom of the plot region to leve room for annotation
; and labels.
      ip=ipanel-1
      if dir eq 0 then begin
         ipx=1 + ip mod ncol
         ipy=nrow - ip/ncol
      endif else begin
         ipx=1 + ip/nrow
         ipy=nrow - ip mod nrow
      endelse
      if frx gt 0.9 then xx=frx-0.9 else xx=0.0
      xw=(1.0-xx)/ncol
      x0=xx+xw*(ipx-frx)
      x1=xx+xw*(ipx)
      if fry gt 0.9 then yy=fry-0.9 else yy=0.0
      yw=(1.0-yy)/nrow
      y0=yy+yw*(ipy-fry)
      y1=yy+yw*(ipy)
      return
      end
