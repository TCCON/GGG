      subroutine subtract_channel_fringe(mmp,nmp,ncbf,resid,calcul,
     & obsrvd)
c  Does a Fourier transform of RESID.
c  Determines the amplitude, frequency & phase of largest channel fringe.
c  Subtracts the resulting sine/cosine wave from OBSRVD.
c
c  Inputs:
c     MMP          I*4  Declared dimension of RESID vector
c     NMP          I*4  Used dimension of vectors
c     NCBF         I*4  Number of Continuum Basis Functions
c     RESID(NMP)   R*4  Residual spectrum
c     CALCUL(NMP)  R*4  Calculted spectrum
c     OBSRVD(NMP)  R*4  Observed spectrum
c
c  Outputs:
c     OBSRVD(NMP)  R*4  Observed spectrum (with CF removed)

      implicit none

      integer*4
     & mmp,nmp,ncbf

      real*4
     & obsrvd(nmp),calcul(nmp),resid(nmp),
     & cfamp,cffreq,cfperiod,cfphase,tcbar,
     & tiny,spi,unity

      parameter(tiny=1.0E-18,spi=3.14159265,unity=1.0)    !do not use the value in const_params.f

      call vdot(calcul,1,unity,0,tcbar,nmp)
      tcbar=tiny+tcbar/nmp
      call fit_channel_fringe(mmp,nmp,ncbf,resid,cfamp,cfperiod,cfphase)
c      write(*,*) 'cf amp, period, phase = ',cfamp,cfperiod,cfphase,nmp
      if(cfamp*cfperiod.gt.0.0) then   ! GCT 2016-05-12
         cffreq=2*spi*float(nmp-1)/cfperiod
         call vramp(resid,1,nmp)
         call vsma(resid,1,cffreq,cfphase,0,resid,1,nmp)
         call vcos(resid,1,resid,1,nmp)
         call vmul(resid,1,calcul,1,resid,1,nmp)
         call vsma(resid,1,-(cfamp/tcbar),obsrvd,1,obsrvd,1,nmp)
      end if
      return
      end

