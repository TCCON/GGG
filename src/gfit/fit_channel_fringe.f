      subroutine fit_channel_fringe(mmp,nmp,ncbf,ybuf,
     & cfamp,cfperiod,cfphase)
c  Solves the equation YBUF(i) = CFAMP.COS[2.Pi.i/CFPERIOD+CFPHASE] 
c  by FFT methods. Then parameterizes the channel fringes present
c  in vector YBUF.
c
c  On Input:
c      MMP       I*4  Declared dimension of YBUF
c      NMP       I*4  Number of elements of YBUF to be analyzed
c      NCBF      I*4  Number of Continuum Basis Finctions
c      YBUF()    R*4  Vector of residuals to be analyzed for channel fringes
c
c  On Output:
c      CFAMP     R*4  Amplitude of dominant cosine (same units as YBUF)
c      CFPERIOD  R*4  Period of dominant cosine (in grid points)
c      CFPHASE   R*4  Phase of dominant cosine (radians)
c      YBUF()    R*4  contains the pure cosine wave
c
c  Notes:
c      1) The data originally in YBUF will be destroyed
c      2) MMP must exceed or equal the next power of 2 above NMP
c      3) Returns CFAMP, CFPHASE, CFFREQ all = 0.0 if fringe
c         characteristics could not be unambiguously determined.
c      4) NCBF is used to prevent competition, between the
c        continuum fitting and channel fringes, to model low
c        frequency undulations in the residuals.
c
      implicit none

      integer*4 mmp,nmp,ncbf,i,imax,k,istart,nfft
      real*4 ybuf(mmp),cfamp,cfperiod,cfphase,
     * tiny,aa,spi,
     & xi,am,dm,amax,ap,dp,phasm,phasi,phasp,
     & denom
      parameter (tiny=1.e-18, spi=3.14159265)    ! do not use local params file !
      data k/0/

      k=k+1
      cfamp=0.0
      cfperiod=0.0
      cfphase=0.0
c
c  Find the smallest power of 2 to accomodate NMP
      nfft=1
      do while (nfft.lt.nmp)
         nfft=nfft+nfft
      end do
      if (nfft.gt.mmp) stop 'fringes.f: Increase parameter MMP'
c      do i=1,nmp
c      write(55,*) k,i,ybuf(i)
c      end do
c
c  Re-arrange YBUF vector to FFT-packed format
      call fft_pack(ybuf,nmp,nfft)
c
c  Perform fast Fourier transform.
      call ffak(ybuf,nfft)
c
c  Search for point IMAX having the peak amplitude in power spectrum.
c  Start at point ISTART, skipping low frequencies that the continuum
c  fitting will deal with anyway.
c      istart=7  ! Will find 3.5 or more periods per interval
c      istart=6  ! Will find 3.0 or more periods per interval
c      istart=5  ! Will find 2.5 or more periods per interval
c      istart=4  ! Will find 2.0 or more periods per interval
      istart=max(3,1+ncbf) 
      imax=0
      amax=0.0
      am=cabs(cmplx(ybuf(2*istart-5),ybuf(2*istart-4)))
      aa=cabs(cmplx(ybuf(2*istart-3),ybuf(2*istart-2)))
      do i=istart,(nfft/2)*3/4   ! 75% of Nyquist
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
      if(abs(phasi-phasm).gt. spi) phasm=phasm+sign(2*spi,phasi-phasm)
      if(abs(phasi-phasp).gt. spi) phasp=phasp+sign(2*spi,phasi-phasp)
      dp=ap
      dm=am
      if( abs(2*phasi-phasp-phasm) .gt. spi/2 ) then
         if( ap*abs(phasi-phasm) .gt. am*abs(phasi-phasp)) then
            phasm=atan2(ybuf(2*(imax-1))/am,-(ybuf(2*(imax-1)-1)/am))
            dm=-(4*am)
         else
            phasp=atan2(ybuf(2*(imax+1))/ap,-(ybuf(2*(imax+1)-1)/ap))
            dp=-(4*ap)
         endif
      endif
      if(abs(phasi-phasm).gt. spi) phasm=phasm+sign(2*spi,phasi-phasm)
      if(abs(phasi-phasp).gt. spi) phasp=phasp+sign(2*spi,phasi-phasp)
c
c  Determine frequency of predominant channel fringe.
      denom=2*amax-dp-dm
      if(abs(denom).le.0.0) return
      xi=0.5*(dp-dm)/denom
      cfperiod=float(nfft)/(xi+imax-1)
c
c  Interpolate peak amp and phase to channel fringe frequency.
      cfphase=phasi+xi*((0.5+xi)*(phasp-phasi)+(0.5-xi)*(phasi-phasm))
c      cfphase=(phasm*am+phasi*amax+phasp*ap)/(am+amax+ap) ! old method
      cfamp=(2*amax+xi*(dp-dm))/nmp
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
