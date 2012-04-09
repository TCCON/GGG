      subroutine write_spt(lun_spt,winfo,sptpath,obsrvd,calcul,
     & cx,ex,startmp,dopp,gint,overcol,pars,sza,obalt,zmin,
     & wrms,frac,pd,ssnmp,nmp,nmpfp,ntg,nfp)

      implicit none
      integer*4 nmp, nmpfp, k, jtg, ntg, nfp, lun_spt, fbc
      real*8 startmp,dopp,gint,freq,sza,obalt,wlimit,frac
      real*4 obsrvd(nmp), calcul(nmp), cntuum, zmin, wrms,rr,ss,
     & cx(nfp), ex(nfp), overcol(ntg), pd(nmpfp,ntg+2),ssnmp(nmp)
      character sptpath*128,pars(ntg)*9,winfo*(*)
c
c  Write .spt file
      open(lun_spt,file=sptpath(:fbc(sptpath)-1),status='unknown',
     & err=66) 
      if( index(winfo,' so ').gt.0 .or. index(winfo,' so/').gt.0) then
        write(lun_spt,*)3,ntg+5  ! Freq, Meas, Calc, Targets, Other, Solar
        write(lun_spt,'(2f14.6,i7,3f8.3,1x,2f7.4,21(1pe11.3,e8.1))',
     &  err=67)
     &  startmp*(1.0d0+dopp),(startmp+gint*(nmp-1))*(1.0d0+dopp),
     &  nmp,sza,obalt,zmin,wlimit(dble(wrms),'f7.4'),frac,
     &  (cx(jtg)*overcol(jtg),ex(jtg)*overcol(jtg),jtg=1,ntg)
        write(lun_spt,'(a40,19a13)')
     &   '   Freq         Tm           Tc         ',
     &   ('    '//pars(jtg),jtg=1,ntg),'   other     ','   solar     '
         do k=1,nmp
            freq=(startmp+(k-1)*gint)*(1.0d0+dopp)
            rr=-0.5+float(k-1)/(nmp-1)
            ss = (3*(nmp-1)*rr**2-0.25*(nmp+1))/(nmp-2)
            cntuum=cx(ntg+1)*(1.+cx(ntg+2)*rr+cx(ntg+3)*ss)
c            cntuum=cx(ntg+1)*(1.+cx(ntg+2)*float(k-(nmp+1)/2)/(nmp-1))
            write(lun_spt,'(f12.6,1p20e13.5)',err=66) freq,
     &      obsrvd(k)/cntuum-cx(ntg+5), calcul(k)/cntuum-cx(ntg+5),
     &      (pd(k,jtg),jtg=2,ntg+1),pd(k,1),ssnmp(k)
         end do
         close(lun_spt)
      else
        write(lun_spt,*)3,ntg+4  ! Freq, Meas, Calc, Targets, Other
        write(lun_spt,'(2f14.6,i7,3f8.3,1x,2f7.4,21(1pe11.3,e8.1))',
     &  err=67)
     &  startmp*(1.0d0+dopp),(startmp+gint*(nmp-1))*(1.0d0+dopp),
     &  nmp,sza,obalt,zmin,wlimit(dble(wrms),'f7.4'),frac,
     &  (cx(jtg)*overcol(jtg),ex(jtg)*overcol(jtg),jtg=1,ntg)
        write(lun_spt,'(a40,18a13)')
     &    '   Freq         Tm           Tc         ',
     &   ('    '//pars(jtg),jtg=1,ntg),'   other     '
         do k=1,nmp
            freq=(startmp+(k-1)*gint)*(1.0d0+dopp)
            cntuum=cx(ntg+1)*(1.+cx(ntg+2)*float(k-(nmp+1)/2)/(nmp-1))
            write(lun_spt,'(f12.6,1p20e13.5)',err=66) freq,
     &      obsrvd(k)/cntuum-cx(ntg+5), calcul(k)/cntuum-cx(ntg+5),
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
