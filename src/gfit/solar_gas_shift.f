      function solar_gas_shift(cl,ct,obsrvd,calcul,spts,nmp)
c
c  code developed 20050113 to calculate the difference between
c  solar and gas shifts.  Modified 20060106.
c  Converted to a subroutine 20070818
c
c  Inputs:
c      CL           R*4  Continuum Level
c      CT           R*4  Continuum Tilt
c      OBSRVD(NMP)  R*4  Measured spectrum
c      CALCUL(NMP)  R*4  Calculated spectrum
c      STS(NMP)     R*4  Solar Pseudo Transmittance Spectrum (at the
c                        grid spacing of the measured spectrum)
c      NMP          I*4  Number of Measured Points
c

      implicit none
      integer*4 k,nmp,nop
      real*4 
     & cl,           ! Continuum level
     & ct,           ! Continuum Tilt
     & obsrvd(nmp),  ! Observed spectrum
     & calcul(nmp),  !  calculated transmittance spectrum
     & spts(nmp),     !  Solar transmittance spectrum
     & res, continuum,
     & gs,gg,tns,tds,tng,tdg,
     & solar_gas_shift ! in units of the measured spectral point spacing

       tns=0.0
       tds=1.e-16
       tng=0.0
       tdg=1.e-16
       nop=(nmp+1)/2
       do k=2,nmp-1
          continuum=cl*(1.+ct*float(k-nop)/(nmp-1))
          res=obsrvd(k)-calcul(k)
c
c Solar Shift: (spts is the solar transmittance spectrum)
          gs=continuum*(spts(k+1)-spts(k-1)) ! partial Differential
          tns=tns+gs*res
          tds=tds+gs**2
c
c Gas shift (calcul/spts is the gaseous transmittance spectrum)
          gg=continuum*(calcul(k+1)/spts(k+1)-calcul(k-1)/spts(k-1))
          tng=tng+gg*res
          tdg=tdg+gg**2
       end do     ! k=2,nmp-1
       solar_gas_shift=tns/tds-tng/tdg
       return
       end
