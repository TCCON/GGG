      subroutine subtract_cf(obsrvd,calcul,resids,nmp,mmp)
c  Does a Fourier transfor of RESIDS
c  Determines the amplitude, frequency & phase of any channel fringes.
c  Subtracts the resulting sine/cosine wave from OBSRVD.

      implicit none
      include "../ggg_const_params.f"

      integer*4
     & mmp,nmp

      real*4
     & obsrvd(nmp),calcul(nmp),resids(nmp),
     & cfamp,cffreq,cfphase,tcbar,
     & tiny

      parameter(tiny=1.0E-18)    !do not use the value in const_params.f
c
      call vdot(calcul,1,unity,0,tcbar,nmp)
      tcbar=tiny+tcbar/nmp
      call fringes(cfamp,cffreq,cfphase,resids,nmp,mmp)
      call vramp(resids,1,nmp)
      call vsma(resids,1,cffreq,cfphase,0,resids,1,nmp)
      call vcos(resids,1,resids,1,nmp)
      call vmul(resids,1,calcul,1,resids,1,nmp)
      call vsma(resids,1,-(cfamp/tcbar),obsrvd,1,obsrvd,1,nmp)
      return
      end

