      subroutine write_spt(lunv,winfo,sptpath,obsrvd,calcul,
     & cx,ex,startmp,dopp,gint,overcol,pars,sza,obalt,zmin,
     & wrms,pd,nmp,nmpfp,ntg)

      integer*4 nmp, nmpfp, k, jtg, ntg, lunv, fbc
      real*8 startmp,dopp, gint, freq, sza, obalt
      real*4 obsrvd(nmp), calcul(nmp), cntuum, zmin, wrms,
     & cx(ntg+4), ex(ntg+4), overcol(ntg), pd(nmpfp,ntg+2)
      character sptpath*128,pars(ntg)*9,winfo*(*)
c
c  Write .spt file
      open(lunv,file=sptpath(:fbc(sptpath)-1),status='unknown',err=66) 
      if( index(winfo,' so ') .gt. 0) then
         write(lunv,*)3,ntg+5
         write(lunv,'(2f14.6,i7,f8.3,2f8.3,f7.4,1pe11.3,e8.1)',err=66)
     &   startmp*(1.0d0+dopp),(startmp+gint*(nmp-1))*(1.0d0+dopp),
     &   nmp,sza,obalt,zmin,wrms,cx(1)*overcol(1),ex(1)*overcol(1)
         write(lunv,'(a40,16a13)')
     &   '   Freq         Tm           Tc         ',
     &   ('    '//pars(jtg),jtg=1,ntg),'   other     ','   solar     '
         do k=1,nmp
            freq=(startmp+(k-1)*gint)*(1.0d0+dopp)
            cntuum=cx(ntg+1)*(1.+cx(ntg+2)*float(k-(nmp+1)/2)/(nmp-1))
            write(lunv,'(f12.6,1p12e13.5)',err=66) freq,
     &      obsrvd(k)/cntuum-cx(ntg+4), calcul(k)/cntuum-cx(ntg+4),
     &      (pd(k,jtg),jtg=2,ntg+1),pd(k,1),pd(k,ntg+2)
         end do
         close(lunv)
      else
         write(lunv,*)3,ntg+4
         write(lunv,'(2f14.6,i7,f8.3,2f8.3,f7.4,1pe11.3,e8.1)',err=66)
     &   startmp*(1.0d0+dopp),(startmp+gint*(nmp-1))*(1.0d0+dopp),
     &   nmp,sza,obalt,zmin,wrms,cx(1)*overcol(1),ex(1)*overcol(1)
         write(lunv,'(a40,16a13)')
     &    '   Freq         Tm           Tc         ',
     &   ('    '//pars(jtg),jtg=1,ntg),'   other     '
         do k=1,nmp
            freq=(startmp+(k-1)*gint)*(1.0d0+dopp)
            cntuum=cx(ntg+1)*(1.+cx(ntg+2)*float(k-(nmp+1)/2)/(nmp-1))
            write(lunv,'(f12.6,1p12e13.5)',err=66) freq,
     &      obsrvd(k)/cntuum-cx(ntg+4), calcul(k)/cntuum-cx(ntg+4),
     &      (pd(k,jtg),jtg=2,ntg+1),pd(k,1)
         end do
      endif
      close(lunv)
      return  ! Normal return

 66   continue
      write (*,*)'Warning: SPT files will not be written'
      write (*,*)sptpath//' does not exist'
      close(lunv)
      return  ! Error return
      end
