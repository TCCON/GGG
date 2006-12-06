      subroutine fringes(cfamp,cffreq,cfphase,resids,nmp,mmp)
c  Solves the equation RESID = CFAMP.COS[CFFREQ.RAMP+CFPHASE] by FFT methods
c  in order to parameterizes the channel fringes present in vector RESID.
c
c  On Input:
c      RESID     R*4  array of residuals to be analyzed for channel fringes
c      NMP       I*4  Number of elements of RESID to be analyzed
c      MMP       I*4  Declared dimension of RESID
c
c  On Output:
c      CFAMP     R*4  Amplitude of dominant cosine
c      CFFREQ    R*4  Frequency of dominant cosine
c      CFPHASE   R*4  Phase of dominant cosine
c      RESID     R*4  contains the pure cosing wave
c
c  Notes:
c      1) The data originally in RESID will be destroyed
c      2) MMP must exceed or equal the next power of 2 above NMP
c      3) Returns CFAMP, CFPHASE, CFFREQ all = 0.0 if fringe
c         characteristics could not be unambiguously determined.
c
      implicit none
      integer*4 nfft,nmp,mmp,nop,i,imax
      real*4 resids(mmp),cfamp,cffreq,cfphase,zero,tiny,aa,
     & xi,am,dm,amax,ap,dp,phasm,phasi,phasp,dd,pi,denom
      parameter (zero=0.0,tiny=0.e-18,pi=3.14159265)
c
      cfamp=0.0
      cffreq=0.0
      cfphase=0.0
c
c  Find the smallest power of 2 to accomodate NMP
      nfft=1
      do while (nfft.lt.nmp)
      nfft=nfft+nfft
      end do
      if (nfft.gt.mmp) stop 'fringes: Increase parameter MMP'
c
c  Re-arrange RESID vector to FFT-packed format
      nop=(nmp-1)/2
      dd=resids(1+2*nop)
      call vswap(resids(1),1,resids(1+nop),1,nop)
      call vmov(resids(2*nop),-1,resids(nfft),-1,nop)
      resids(nop+1)=dd
      call vmov(zero,0,resids(nop+2),1,nfft-2*nop-1)
c
c  Perform fast Fourier transform.
      call ffak(resids,nfft)
c
c  Search for point IMAX having the peak amp (ignoring DC and Nyquist terms)
      imax=0
      amax=0.0
      am=abs(resids(1))
      aa=cabs(cmplx(resids(3),resids(4)))
      do i=3,nfft/2
        ap=cabs(cmplx(resids(2*i-1),resids(2*i)))
        if(aa.ge.amax) then
        if(aa.ge.am .and. aa.ge.ap) then
           amax=aa
           imax=i-1
        endif
        endif
        am=aa
        aa=ap
      end do
c
c  Test to see if a true maximum was found.
      if(imax.eq.0) return
c
c  Use point IMAX, together with its two neighbors, to determine the
c  fractional frequency point having the peak amp. 
      amax=2*tiny+amax
      am=tiny+cabs(cmplx(resids(2*(imax-1)-1),resids(2*(imax-1))))
      ap=tiny+cabs(cmplx(resids(2*(imax+1)-1),resids(2*(imax+1))))
      phasm=atan2(-resids(2*(imax-1))/am,resids(2*(imax-1)-1)/am)
      phasi=atan2(-resids(2*(imax+0))/amax,resids(2*(imax+0)-1)/amax)
      phasp=atan2(-resids(2*(imax+1))/ap,resids(2*(imax+1)-1)/ap)
      if(abs(phasi-phasm).gt. pi) phasm=phasm+sign(2*pi,phasi-phasm)
      if(abs(phasi-phasp).gt. pi) phasp=phasp+sign(2*pi,phasi-phasp)
      dp=ap
      dm=am
      if( abs(2*phasi-phasp-phasm) .gt. pi/2 ) then
      if( ap*abs(phasi-phasm) .gt. am*abs(phasi-phasp)) then
        phasm=atan2(resids(2*(imax-1))/am,-resids(2*(imax-1)-1)/am)
        dm=-4*am
      else
        phasp=atan2(resids(2*(imax+1))/ap,-resids(2*(imax+1)-1)/ap)
        dp=-4*ap
      endif
      endif
      if(abs(phasi-phasm).gt. pi) phasm=phasm+sign(2*pi,phasi-phasm)
      if(abs(phasi-phasp).gt. pi) phasp=phasp+sign(2*pi,phasi-phasp)
c
c  Determine frequency of predominant channel fringe.
      denom=2*amax-dp-dm
      if(denom.eq.0.0) return
      xi=0.5*(dp-dm)/denom
      cffreq=2*pi*float(nmp-1)*(xi+imax-1)/nfft
c
c  Interpolate peak amp and phase to channel fringe frequency.
      cfphase=phasi+xi*((0.5+xi)*(phasp-phasi)+(0.5-xi)*(phasi-phasm))
c      cfphase=(phasm*am+phasi*amax+phasp*ap)/(am+amax+ap) ! old method
      cfamp=(2*amax+xi*(dp-dm))/nmp
c      write(*,'(4f12.4)')xi,cfamp,cffreq,cfphase
      return
      end
