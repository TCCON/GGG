      subroutine write_spt(lun_spt,flagso,sptpath,obsrvd,calcul,
     & cont,cx,ex,startmp,dopp,effres,gint,overcol,pars,sza,obalt,zmin,
     & xzo,rmsocl,peff,frac,pd,ssnmp,nmp,ntg,nfp)

c  Inputs:
c     lun_spt      I*4   Logical Unit Number 
c     flagso       I*4   1 means solar spectrum, 0 otherwise
c     sptpath      C*128 Output file path
c     obsrvd(nmp)  R*4   Observed spectrum (Arbitrary units)
c     calcul(nmp)  R*4   Calculated spectrum (Arbitrary units)
c                        Calcul = Cont*(SPTS*TRAN*(1-XZO)+XZO)
c     cont(nmp)    R*4   Continuum spectrum
c     cx(nfp)      R*4   State-Vector: current values
c     ex(nfp)      R*4   State-Vector: uncertainties
c     startmp      R(8   wavenumber of first point
c     dopp         R*8   doppler stretch
c     effres       R*8   Effective Spectral Resolution (cm-1)
c     gint         R*8   Spectral point spacing (cm-1)
c     overcol(ntg) R*8   vertical columns of targest gases
c     pars(ntg)    C*9   names of targes gases
c     sza          R*8   Solar Zenith angle
c     obalt        R*8   Observation Altitude (km)
c     zmin         R*4   Minimum altitude
c     xzo          R*4   Zero Offset
c     rmsocl       R*4   rms spectral fit
c     peff         R*4   Effective pressure (first target gas)
c     frac         R*4   Fraction of solar disk observed
c     pd(nmp,ntg)  R*4   Individual target gas transmittances
c     ssnmp(nmp)   R*4   Solar pseudo-transmittance spectrum
c     nmp          I*4   Number of Measured (spectral) Points
c     ntg          I*4   Number of target gases
c     nfp          I*4   Number of fitted parameters = ntg+cl+ct+fs+?
c
c  Output:
c     output file written to SPTPATH

      implicit none
      integer*4 nmp, j, k, jtg, ntg, nfp, lun_spt, fbc,flagso
      real*8 startmp,dopp,effres,gint,freq,sza,obalt,wlimit,frac
      real*4 obsrvd(nmp),calcul(nmp),zmin,rmsocl,cont(nmp),xzo,peff,
     & cx(nfp), ex(nfp), overcol(ntg), pd(nmp+nfp,ntg+2),ssnmp(nmp)
      character sptpath*128,pars(ntg)*(*)

      open(lun_spt,file=sptpath(:fbc(sptpath)-1),status='unknown',
     & err=66) 
      write(lun_spt,*)3,ntg+5+flagso  ! Freq, Meas, Calc, Cont, Targets, Other, Solar
      write(lun_spt,
     &'(2f14.6,i7,f8.4,3f8.3,1x,f9.6,1pe10.3,0pf7.4,f8.4,20(1p3e11.3))',
     &  err=67) startmp*(1.d0+dopp),(startmp+gint*(nmp-1))*(1.d0+dopp),
     &  nmp,effres,sza,obalt,zmin,wlimit(dabs(dble(rmsocl)),'f8.6'),
     &  peff,frac,xzo,(overcol(jtg),cx(jtg),ex(jtg),jtg=1,ntg)
      write(lun_spt,'(a50,19a13)')
     &   '   Freq         Tm           Tc           Cont    ',
     &   ('     '//pars(jtg),jtg=1,ntg),'   other     ',
     &   ('   solar',j=1,flagso)
c  SPTS*TRAN = [Calcul/Cont-XZO]/(1-XZO)
      do k=1,nmp
         freq=(startmp+(k-1)*gint)*(1.0d0+dopp)
         write(lun_spt,'(f12.6,1p21e13.5)',err=66) freq,
c     &   (obsrvd(k)/cont(k)-xzo)/(1-xzo),
c     &   (calcul(k)/cont(k)-xzo)/(1-xzo),cont(k),
     &   obsrvd(k),calcul(k),cont(k),
     &   (pd(k,jtg),jtg=2,ntg+1),pd(k,1),(ssnmp(k),j=1,flagso)
      end do
      close(lun_spt)
      return  ! Normal return

 66   close(lun_spt)
      write (*,*)'write_spt: Error opening .spt file: ',sptpath
      return  ! Error return
 67   close(lun_spt)
      write (*,*)'write_spt: Error writing .spt file: ',sptpath
      return  ! Error return
      end
