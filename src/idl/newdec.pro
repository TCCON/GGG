pro newdec,fin,nin,oper,nop,intrp,fout,nout
;
; Test for cases which will cause array-bound violations
      nhw=(nop-1)/intrp/2
      if nhw lt 1  then  print,' NEWDEC called with NOP < 1+2*ODEC'
      khi=(nout-1)/intrp+1+2*nhw
      if khi gt nin then   print,' NEWDEC warning: KHI > NIN:',khi,nin
;
;  Main loop
      ix=long(intrp)
      for kout=long(0),long(nout-1) do begin
         nn=ix-1
         kin=1+nn/intrp
         kop=intrp - (nn mod intrp)
;         call vdot(fin(kin),1,oper(kop+0),intrp,fout(kout),2*nhw)
         dp=0.0
         for j=0,2*nhw-1 do begin
             dp=dp+fin(kin-1)*oper(kop-1)
             kin=kin+1
             kop=kop+intrp
         endfor
         fout(kout)=dp
         ix=ix+1
      endfor    ;  do kout=0,nout-1
      return
      end
