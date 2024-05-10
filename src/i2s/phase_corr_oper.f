      subroutine phase_corr_oper(lpco,pco,dfr,pco_thresh,phasepath)
c  Takes the interferogram center burst (PCO) and calculates an operator that,
c  when convolved with the raw interferogram, should make it symmetrical.
c 
c  First apodizes the raw interferogram center burst in array PCO(LPCO).
c  Then computes the complex low resolution spectrum (a+i.b) where i=sqrt(-1).
c  The complex low-resolution spectrum contains LPCO/2+1 real values (a) and
c  LPCO/2-1 imaginary values (b). The first and last b values, representing
c  the DC and Nyquist frequencies, are zero by definition because the DC term
c  is a constant value, which is symmetrical about any point. And a Nyquist-
c  sampled function is +/- a constant value, which is symmetrical about any point.
c  These values are omitted from the returned spectrum and the first imaginary
c  element of the returned spectrum actually contains the real Nyquist amplitude,
c  not the DC imaginary.  This is commonly called "FFT packing".
c
c  Computes phase = arctan2(b/a), interpolating across regions with low power.
c  Then computes the complex exponential of the adjusted phase spectrum,
c  and performs a complex-to-real FFT to yield the phase correction operator.
c
c  Note, if the original phase were unadjusted, its complex exponential would
c  equal  (a-i.b)/sqrt(a^2+b^2). Therefore, multiplying the complex spectrum by 
c  it results in a spectrum = (sqrt(a^2+b^2),0) with imaginary terms all zero.
c  Thus, the complex-to-real FFT of this spectral product, the phase corrected
c  interferogram, would be perfectly symmetrical.
c
c  INPUTS:
c     LPCO           Length of Phase Correction Operator
c     PCO(LPCO)      Array containing FFT-packed interferogram center burst
c     DFR            Fractional displacement of exact (non-integer) ZPD
c     pco_thresh     Fractional amplitude threshold for phase information
c     phasepath      Where to write the phase curve
c
c OUTPUTS:
c     PCO(LPCO)       Array containing FFT-packed phase correction operator

      implicit none

      integer*4 mpco,lpco,i,lastoki,lnbc
      parameter (mpco=64*1024)
      real*4 pco(lpco),amplit(1+mpco/2),phase(1+mpco/2),area,rarea,
     & ampmax,pythag
c      real*4 pco_odd(mpco/2),pco_even(mpco/2)
      real*8 pco_thresh,flag,freq,dfr
      character phasepath*(128)
c
      if(lpco.gt.mpco) stop 'Error in phase_corr_oper: LPCO > MPCO'

      lastoki=0
      flag=1.0d0
      if(pco(1).lt.0.0) flag=-1.d0  ! If ZPD igram value -ve

c  Apply cos^2 windowing igram center-burst in PCO
c  in order to apodize resulting low-res spectrum.
      call window(lpco,pco,dfr)

c  Perform small real-to-complex FFT of igram center burst
      call ffak ( pco, lpco )

c  Compute low resolution power and phase spectra,
c      tot=0.0d0
c      toty=0.0d0
      pco(1)=abs(pco(1))     ! Set DC magnitude 
      pco(2)=0.0            ! Set DC phase 
      amplit(1)=pco(1)      ! DC
      ampmax=amplit(1)
      ampmax=0.0
      phase(1)=0.0         ! DC
      do i=2,lpco/2
         freq=15798*float(i-1)/(lpco/2)
         amplit(i)=pythag(pco(2*i-1),pco(2*i))
         if(amplit(i).gt.ampmax) ampmax=amplit(i)
         if(amplit(i).gt.pco_thresh*ampmax) lastoki=i
         phase(i)=sngl( datan2 ((-flag)*pco(2*i), flag*pco(2*i-1) ))
         pco(2*i-1)=amplit(i)
         pco(2*i)=phase(i)
      end do
c
c  Interpolate phase across regions with little power (mag < pco_thresh*ampmax)
      if(lastoki.eq.0) then
         write(*,*)'Warning: lastoki being set to lpco/2. This is '//
     &  'probably a bad spectrum.'
         lastoki=lpco/2
      endif
      call vintrp(phase,1,lpco/2,lastoki,1)  ! Extrapolate phase to 0 at high wav
      lastoki=1
      do i=lpco/32,lpco/2
         if(amplit(i).gt.pco_thresh*ampmax) then
            if(i.ne.lastoki+1) call vintrp(phase,1,lpco/2,lastoki,i)
            lastoki=i
         endif
      end do
c
c  Write out phase and amplitude curves
      if(lnbc(phasepath).gt.0) then
         open(19,file=phasepath)
         write(19,*)2,6
         write(19,*)' Wavenumber Amplitude Log[Amplitude] Log[Thresh]
     &   Phase1 Phase2'
         do i=1,lpco/2
            freq=15798*float(i-1)/(lpco/2)
            write(19,'(f9.1,5f9.4)')freq,amplit(i),log10(amplit(i)+1E-6)
     &      ,log10(pco_thresh*ampmax+1e-6),pco(2*i),phase(i)
         end do
         close(19)
      endif
c
c  Compute complex exponential of the adjusted phase
      do i=1,lpco/2
         pco(2*i-1)=cos(phase(i))
         pco(2*i)=sin(phase(i))
      end do
      pco(2)=1.0   ! insert the Nyquist frequency
c
c  Take the inverse transform (complex-to-real)
      call ffsk ( pco, lpco )
c
c  Window phase correction operator
      call window( lpco, pco, dfr )
c
c  Normalize the phase-correcting function to unit area
      call vdot(pco,1,1.0,0,area,lpco)
      rarea=sngl(flag)/abs(area)
      call vmul(pco,1,rarea,0,pco,1,lpco)
      return
      end

      subroutine window(lpco,y,dfr)
c   Performs in-place windowing of the FFT-packed igram in
c   vector Y  using the function cos((i-1-fr)*pi/LPCO)^2
c   This will reduce ringing in spectrum produced by FFT of Y.
c
c   Assumes that igram Y is already FFT-packed, i.e. with the
c   +ve half of operator in  Y(1:LPCO/2) and
c   -ve half in              Y(1+LPCO/2:LPCO)
c   such that the nearest point to ZPD is at Y(1).
c
c  Inputs:
c      LPCO    I*4    Length of Phase Correction Operator
c     Y(LPCO)  R*4    Interferogram center-burst
c       FR     R*8    Displacement of exact ZPD location from 1

      implicit none

      integer*4 lpco,i
      real*4 y(lpco)
      real*8 pi,q,x,dfr
c
      parameter (pi = 4.d0*datan(1.d0))

      q=pi/lpco
      x=-q*dfr
      do i=1,lpco
c         write(*,*)i,x,sngl(dcos(x)**2),y(i)
         y(i)=y(i)*sngl(dcos(x)**2)
         x=x+q
      enddo
      return
      end

      subroutine vintrp(y,incr,nn,i1,i2)
      implicit none
      integer incr,nn,i1,i2,j,k,l
      real y(incr,nn),step
      k=i2-i1
      if(k.le.0) k=k+nn
      step=(y(incr,i2)-y(incr,i1))/float(k)
      do  j=1,k-1
         l=mod(i1+j-1,nn)+1
         y(incr,l)=y(incr,i1)+float(j)*step
      enddo
      return
      end
