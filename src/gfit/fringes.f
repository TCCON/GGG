      subroutine fringes(cfamp,cffreq,cfphase,ybuf,nmp,mmp)
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
c      RESID     R*4  contains the pure cosine wave
c
c  Notes:
c      1) The data originally in RESID will be destroyed
c      2) MMP must exceed or equal the next power of 2 above NMP
c      3) Returns CFAMP, CFPHASE, CFFREQ all = 0.0 if fringe
c         characteristics could not be unambiguously determined.
c
      implicit none
      integer*4 nfft,nmp,mmp,i,imax
      real*4 ybuf(mmp),cfamp,cffreq,cfphase,zero,tiny,aa,
     & xi,am,dm,amax,ap,dp,phasm,phasi,phasp,pi,denom
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
      call fft_pack(ybuf,nmp,nfft)
c      nop=(nmp-1)/2
c      dd=ybuf(1+2*nop)
c      call vswap(ybuf(1),1,ybuf(1+nop),1,nop)
c      call vmov(ybuf(2*nop),-1,ybuf(nfft),-1,nop)
c      ybuf(nop+1)=dd
c      call vmov(zero,0,ybuf(nop+2),1,nfft-2*nop-1)
c
c  Perform fast Fourier transform.
      call ffak(ybuf,nfft)
c
c  Search for point IMAX having the peak amp (ignoring DC and Nyquist terms)
      imax=0
      amax=0.0
      am=abs(ybuf(1))
      aa=cabs(cmplx(ybuf(3),ybuf(4)))
c      do i=3,nfft/2
      do i=3,(nfft/2) * 3/4   ! 75% of Nyquist  GCT 2010-07-22
        ap=cabs(cmplx(ybuf(2*i-1),ybuf(2*i)))
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
      am=tiny+cabs(cmplx(ybuf(2*(imax-1)-1),ybuf(2*(imax-1))))
      ap=tiny+cabs(cmplx(ybuf(2*(imax+1)-1),ybuf(2*(imax+1))))
      phasm=atan2(-(ybuf(2*(imax-1))/am),ybuf(2*(imax-1)-1)/am)
      phasi=atan2(-(ybuf(2*(imax+0))/amax),ybuf(2*(imax+0)-1)/amax)
      phasp=atan2(-(ybuf(2*(imax+1))/ap),ybuf(2*(imax+1)-1)/ap)
      if(abs(phasi-phasm).gt. pi) phasm=phasm+sign(2*pi,phasi-phasm)
      if(abs(phasi-phasp).gt. pi) phasp=phasp+sign(2*pi,phasi-phasp)
      dp=ap
      dm=am
      if( abs(2*phasi-phasp-phasm) .gt. pi/2 ) then
      if( ap*abs(phasi-phasm) .gt. am*abs(phasi-phasp)) then
        phasm=atan2(ybuf(2*(imax-1))/am,-(ybuf(2*(imax-1)-1)/am))
        dm=-(4*am)
      else
        phasp=atan2(ybuf(2*(imax+1))/ap,-(ybuf(2*(imax+1)-1)/ap))
        dp=-(4*ap)
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
c      write(*,'(a8,2i6,4f12.4)')'fringes:',nfft,imax,xi,
c     & cfamp,cffreq,cfphase
      return
      end


      subroutine fft_pack(ybuf,nmp,nfft)
c  Performs an in-place FFT-packing of vector Y
c  This involves moving the samples from
c  NMP/2-2 to NFFT-2
c  NMP/2-1 to NFFT-1
c  NMP/2+0 to NFFT
c  NMP/2+1 to 1 
c  NMP/2+2 to 2 
c  NMP/2+3 to 3 
c  etc, etc,
c  It fills in the middle NFFT-NMP points with zeros.
c
c  For the special case where nmp=nfft, fft-pack simply
c  swaps the first and second halves of the y-vector.
c
      integer*4 nmp,nfft,nop,i,k
      real*4 ybuf(nfft),temp,ylast

      nop=nmp/2
      k=nfft
      ylast=ybuf(nmp)
      do i=nop,1,-1
         temp=ybuf(i)
         ybuf(i)=ybuf(i+nop)
         ybuf(k)=temp
         k=k-1
      end do
      do i=1+nop,nfft-nop
         ybuf(i)=0.0
      end do
      if(2*nop.lt.nmp) ybuf(nop+1)=ylast    ! NMP odd
      return
      end
