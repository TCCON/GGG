      subroutine write_spt(lun_spt,winfo,sptpath,obsrvd,calcul,
     & cont,cx,ex,startmp,dopp,gint,overcol,pars,sza,obalt,zmin,
     & xzo,wrms,frac,pd,ssnmp,nmp,nmpfp,ntg,nfp)

      implicit none
      integer*4 nmp, nmpfp, k, jtg, ntg, nfp, lun_spt, fbc
      real*8 startmp,dopp,gint,freq,sza,obalt,wlimit,frac
      real*4 obsrvd(nmp), calcul(nmp), zmin, wrms, cont(nmp), xzo,
     & cx(nfp), ex(nfp), overcol(ntg), pd(nmpfp,ntg+2),ssnmp(nmp)
      character sptpath*128,pars(ntg)*9,winfo*(*)
c
c  Write .spt file
c      if(ipzo.gt.0) then
c         xzo=cx(ipzo)
c      else
c         xzo=0.0
c      endif
      open(lun_spt,file=sptpath(:fbc(sptpath)-1),status='unknown',
     & err=66) 
      if( index(winfo,' so ').gt.0 .or. index(winfo,' so/').gt.0) then
        write(lun_spt,*)3,ntg+5  ! Freq, Meas, Calc, Targets, Other, Solar
        write(lun_spt,'(2f14.6,i7,3f8.3,1x,2f7.4,21(1p3e11.3))',err=67)
     &  startmp*(1.0d0+dopp),(startmp+gint*(nmp-1))*(1.0d0+dopp),
     &  nmp,sza,obalt,zmin,wlimit(dble(wrms),'f7.4'),frac,
     &  (overcol(jtg),cx(jtg),ex(jtg),jtg=1,ntg)
        write(lun_spt,'(a40,19a13)')
     &   '   Freq         Tm           Tc         ',
     &   ('    '//pars(jtg),jtg=1,ntg),'   other     ','   solar     '
         do k=1,nmp
            freq=(startmp+(k-1)*gint)*(1.0d0+dopp)
            write(lun_spt,'(f12.6,1p20e13.5)',err=66) freq,
     &      obsrvd(k)/cont(k)-xzo, calcul(k)/cont(k)-xzo,
     &      (pd(k,jtg),jtg=2,ntg+1),pd(k,1),ssnmp(k)
         end do
         close(lun_spt)
      else
        write(lun_spt,*)3,ntg+4  ! Freq, Meas, Calc, Targets, Other
        write(lun_spt,'(2f14.6,i7,3f8.3,1x,2f7.4,21(1p3e11.3))',err=67)
     &  startmp*(1.0d0+dopp),(startmp+gint*(nmp-1))*(1.0d0+dopp),
     &  nmp,sza,obalt,zmin,wlimit(dble(wrms),'f7.4'),frac,
     &  (overcol(jtg),cx(jtg),ex(jtg),jtg=1,ntg)
        write(lun_spt,'(a40,18a13)')
     &    '   Freq         Tm           Tc         ',
     &   ('    '//pars(jtg),jtg=1,ntg),'   other     '
         do k=1,nmp
            freq=(startmp+(k-1)*gint)*(1.0d0+dopp)
            write(lun_spt,'(f12.6,1p20e13.5)',err=66) freq,
     &      obsrvd(k)/cont(k)-xzo, calcul(k)/cont(k)-xzo,
     &      (pd(k,jtg),jtg=2,ntg+1),pd(k,1)
         end do
      endif
      close(lun_spt)
      return  ! Normal return

 66   close(lun_spt)
      write (*,*)'write_spt: Error opening .spt file',sptpath
      return  ! Error return
 67   close(lun_spt)
      write (*,*)'write_spt: Error writing .spt file',sptpath
      return  ! Error return
      end
