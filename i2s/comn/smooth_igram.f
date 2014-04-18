      subroutine smooth_igram(nip,yin,yout,frsp,dclevel,frzpda,
     & pinv,pinl)
c
c  Smoothes an interferogram by direct convolution with an
c  internally-computed low-pass operator.
c
c  Inputs:
c           nip         I*4      Number of interferogram points
c           yin(nip)    R*4      Raw Interferogram
c           frsp        R*4      Fraction of spectrum passed by filter
c
c Outputs: 
c           yout(nip)   R*4      Smoothed Interferogram
c           dclevel     R*8      Smoothed igram signal level at PINL
c           frzpda      R*8      Fractional amplitude of ZPD dip
c           pinv        R*4      Peak AC interferogram value
c           pinl        I*4      Peak AC interferogram location
c
c Choice of input parameters
c   frsp is the fraction of the frequency domain to be unsmoothed.
c   The lower FRSP, the smoother the resulting interferogram
c   For example, if the spectrum extends from 0 to 15798 cm-1,
c   and you are using an extended InGaAs detector with a cut-on
c   frequency of 3900 cm-1, then set frsp = 3900/15798=0.25
c   For an Si spectrum covering 0 to 15798 cm-1, there is no
c   signal below 9000 cm-1, so set  frsp = 9000/15798=0.57

c   Use the fsrp value that is as large as possible, without
c   intruding into the frequency domain that contains the
c   optical energy of eventual interest. The larger frsp,
c   the narrower the convolution operator that is needed to
c   achieve a given accuracy and therefore the faster the
c   siv_correction process. Choosing the largest frsp-value
c   also allows correction of the most rapid changes in solar
c   intensity. Do not use an artificially small frsp value.
c   This will require a longer convolution operator to achieve
c   a given accuracy.

c   Since the low frequencies that distort the interferogram are
c   generally ~1Hz, whereas the interferometric frequencies are
c   ~1000 Hz, our spectral filter that separates the two can be
c   very gradual/gentle, which means that its Fourier transform
c   (the convolution operator) can be very narrow.
c
c   NHW is the half-width of the convolution operator.
c   An NHW value of 8/frsp achieves a 6-order-of-magnitude
c   suppression of the ZPD. Larger values doesn't really help
c   much in practice because you run into the single precision
c   machine precision.
c
c   The table below shows ZPD supression values achieved by this subroutine
c   for a spectral bandpass extending from 0.35 to 0.65 of fmax
c     NHW=  5      10      15      20      25      30      35      40
c   FRSP ---------------------------------------------------------------
c   0.05 2.6E-03 6.9E-04 2.3E-04 7.8E-05 2.5E-05 6.3E-06 1.4E-06 1.6E-06
c   0.10 1.6E-03 2.0E-04 2.3E-05 1.0E-05 2.5E-06 1.8E-06 6.2E-07 5.5E-07
c   0.15 7.5E-04 8.0E-05 1.5E-05 4.9E-06 2.3E-06 1.1E-06 5.4E-07 3.9E-07
c   0.20 4.3E-04 6.2E-05 1.6E-05 5.2E-06 2.1E-06 1.1E-06 5.7E-07 4.0E-07
c   0.25 6.2E-04 8.3E-05 2.3E-05 8.1E-06 3.4E-06 1.6E-06 8.9E-07 5.4E-07
c   0.30 1.2E-03 2.1E-04 6.3E-05 2.4E-05 1.1E-05 5.3E-06 2.9E-06 1.7E-06
c   0.35 6.5E-03 2.4E-03 1.2E-03 7.7E-04 5.5E-04 4.3E-04 3.5E-04 3.0E-04
c   0.40 1.9E-02 1.1E-02 8.6E-03 7.3E-03 6.6E-03 6.1E-03 5.9E-03 5.7E-03
c   0.45 4.2E-02 3.2E-02 2.8E-02 2.6E-02 2.5E-02 2.5E-02 2.4E-02 2.4E-02
c   0.50 7.3E-02 6.3E-02 6.0E-02 5.8E-02 5.7E-02 5.6E-02 5.6E-02 5.6E-02
c   0.55 1.1E-01 1.0E-01 1.0E-01 9.9E-02 9.8E-02 9.8E-02 9.7E-02 9.7E-02
c   0.60 1.6E-01 1.5E-01 1.5E-01 1.5E-01 1.4E-01 1.4E-01 1.4E-01 1.4E-01
c   0.65 2.1E-01 2.0E-01 2.0E-01 2.0E-01 1.9E-01 1.9E-01 1.9E-01 1.9E-01
c   0.70 2.6E-01 2.5E-01 2.5E-01 2.5E-01 2.5E-01 2.5E-01 2.5E-01 2.5E-01
c   0.75 3.1E-01 3.0E-01 3.0E-01 3.0E-01 3.0E-01 3.0E-01 3.0E-01 3.0E-01
c   0.80 3.5E-01 3.5E-01 3.5E-01 3.5E-01 3.5E-01 3.5E-01 3.5E-01 3.5E-01
c   0.85 4.0E-01 4.0E-01 3.9E-01 3.9E-01 3.9E-01 3.9E-01 3.9E-01 3.9E-01
c   0.90 4.4E-01 4.4E-01 4.4E-01 4.4E-01 4.4E-01 4.4E-01 4.4E-01 4.4E-01
c   0.95 4.8E-01 4.8E-01 4.8E-01 4.8E-01 4.8E-01 4.8E-01 4.8E-01 4.8E-01
c   1.00 5.2E-01 5.2E-01 5.2E-01 5.1E-01 5.1E-01 5.1E-01 5.1E-01 5.1E-01
c-----------------------------------------------------------------------

c    NHW= 4/FRSP  5/FRSP  6/FRSP  7/FRSP  8/FRSP  9/FRSP 10/FRSP 11/FRSP
c   FRSP  --------------------------------------------------------------
c   0.05 2.7E-07 2.7E-07 2.7E-07 2.7E-07 2.7E-07 2.7E-07 2.7E-07 2.7E-07
c   0.10 6.2E-07 3.4E-07 2.8E-07 2.7E-07 2.7E-07 2.7E-07 2.7E-07 2.7E-07
c   0.15 3.3E-06 1.2E-06 5.4E-07 3.5E-07 3.1E-07 2.9E-07 2.8E-07 2.7E-07
c   0.20 1.6E-05 5.2E-06 2.1E-06 1.1E-06 5.7E-07 4.0E-07 3.4E-07 3.6E-07
c   0.25 6.4E-05 2.3E-05 9.7E-06 4.7E-06 2.5E-06 1.4E-06 8.9E-07 6.1E-07
c   0.30 3.8E-04 1.2E-04 6.3E-05 3.4E-05 1.7E-05 1.1E-05 6.9E-06 4.1E-06
c   0.35 5.1E-03 2.8E-03 1.8E-03 1.2E-03 9.1E-04 7.2E-04 5.8E-04 5.2E-04
c   0.40 1.9E-02 1.5E-02 1.1E-02 1.0E-02 8.6E-03 8.0E-03 7.3E-03 7.0E-03
c   0.45 4.6E-02 3.9E-02 3.4E-02 3.0E-02 2.9E-02 2.8E-02 2.7E-02 2.6E-02
c   0.50 8.4E-02 7.3E-02 6.8E-02 6.4E-02 6.2E-02 6.1E-02 6.0E-02 5.9E-02
c   0.55 1.3E-01 1.2E-01 1.1E-01 1.1E-01 1.0E-01 1.0E-01 1.0E-01 1.0E-01
c   0.60 1.7E-01 1.7E-01 1.6E-01 1.5E-01 1.5E-01 1.5E-01 1.5E-01 1.5E-01
c   0.65 2.3E-01 2.1E-01 2.1E-01 2.0E-01 2.0E-01 2.0E-01 2.0E-01 2.0E-01
c   0.70 2.8E-01 2.7E-01 2.6E-01 2.6E-01 2.5E-01 2.5E-01 2.5E-01 2.5E-01
c   0.75 3.3E-01 3.2E-01 3.1E-01 3.1E-01 3.0E-01 3.0E-01 3.0E-01 3.0E-01
c   0.80 3.8E-01 3.7E-01 3.6E-01 3.6E-01 3.5E-01 3.5E-01 3.5E-01 3.5E-01
c   0.85 4.2E-01 4.1E-01 4.1E-01 4.0E-01 4.0E-01 4.0E-01 4.0E-01 4.0E-01
c   0.90 4.7E-01 4.5E-01 4.5E-01 4.4E-01 4.4E-01 4.4E-01 4.4E-01 4.4E-01
c   0.95 5.0E-01 4.9E-01 4.9E-01 4.9E-01 4.8E-01 4.8E-01 4.8E-01 4.8E-01
c   1.00 5.4E-01 5.3E-01 5.2E-01 5.2E-01 5.2E-01 5.2E-01 5.2E-01 5.2E-01
c   The subroutine essentially applies a low-pass filter of the form
c     f(nu) = cos(pi*nu/(2*frsp*dnu))**2     nu < frsp*dnu
c     f(nu) = 0.0                            nu > frsp*dnu
c   to the frequency content of the interferogram (i.e. the spectrum).
c   dnu is the full spectral bandwidth of the raw data.
c
c  In theory, this filter transmits nothing above frsp*dnu
c  and has a 50% transmission at a frequency of 0.5*frsp*dnu
c  In practice, its implementation via a finite convolution
c  operator causes a small leakage of high frequencies due
c  to truncation of the smoothing operator.
c
c  The function f(nu) is known as the Hann apodization function.
c  The smoothing operator is its Fourier transform:
c    F(x') =  sinc(x') + 0.5 * [ sinc(x'-pi) + sinc(x'+pi) ]
c  where x=2*Pi*x*frsp*dnu
c  This is the sum of three sinc functions, the outer
c  ones being half the amplitude of, and centered on the
c  first zeros of, the center sinc function.
c
c  Note that F(x) can be simplified to sinc(x)/[1-(x/pi)**2]
c  We don't actually implement it this way because of the
c  indeterminacies at x = +-pi.  But this form of the equation
c  clearly shows that F(x) falls off as 1/x^3, or 18dB/octave
c  in signal processing terminology, which is much faster than
c  the sinc function which falls of as 1/x or 6 dB/octave.
c
c  So a relatively short convolution operator of the form F(x)
c  can achieve high accuracy. For example, just 20 half-widths
c  (10 full widths) from the center of the operator (x=20.pi),
c  its value has fallen by a factor 25000. So truncating the
c  operator beyond 20 half-widths has only a minor impact, 
c  even less if it is first apodized.
c
c  Convolving the interferogram with F(x) has the same
c  effect as multiplying its Fourier transform by f(nu)
c  and then inverse transforming. For fairly short operator
c  lengths < 2.log2(NFFT) the direct convolution will be
c  faster than the FFT-based convolution.
c
c  When you convolve a NIP-point interferogram with
c  a 2*NHW+1-point operator, the resulting smoothed
c  interferogram will have a length NIP-2*NHW.
c  This is because a smoothing cannot be performed
c  for the first and last NHW points. Although this
c  is not a big deal in terms of information loss,
c  since generally NIP>>NHW, it would be a potential
c  source of confusion if the smoothed interferogram
c  were slightly shorter than the original.
c  The subroutine therefore linearly extrapolates the
c  smoothing function into the "end-zones" to allow
c  the output of a full-size smoothed interferogram.
c  NB: the FFT method of convolution avoids this problem
c  by assuming that the interferogram is a section of an
c  infinite periodic function.
c
c  This subroutine has three advantages over Gretchen's
c  original FFT-based dc-correction subroutine: 
c  1) You don't need a full-size vector in which to store the
c     smoothed interferogram. Only a small workspace array (OPER)
c     is needed
c  2) It is faster for typical operator lengths < 2.log2(NFFT)
c  3) Conceptual simplicity: we're simply smoothing the
c     interferogram by convolving it with an operator.
c     No FFT's required -- no calls to external subroutines
c  
      implicit none
      integer*4
     & nip,nhw,nwid,i,j,pinl,ii

      real*4 
     & yin(nip),       !   Input Data vector
     & yout(nip),      !   Output data vector
     & yac,frsp,pinv,thresh,grad,fwt

      real*8 sinc,xx,pi,smoo,top,ds
      real*8 t0,t1,t2,ty0,ty1,det,wt,dclevel,frzpda

      pi=4*datan(1.0d0)

c  Empirically-determined optimum operator width
c  For large NIP, Min value NHW = 9./1.015 = 9
c  For large NIP, Max value NHW = 9./0.015 = 600
      nhw=min0(nip/2,nint(9./(0.015+abs(frsp)))) ! Max = 9/.015 = 600
c      write(*,*) 'frsp,nhw=',frsp,nhw

c  Pre-compute half of the convolution operatator and weakly apodize,
c  to minimize truncation error. Since the operator is symmetrical,
c  and the center point has a value of 1.0, it is necessary only to
c  save NHW points, which we place in elements 1 thru NHW of YOUT,
c  which otherwise don't get addressed during the convolution.
      top=0.5d0
      do j=1,nhw
         xx=(0.01+abs(frsp))*dfloat(j)*pi
         yout(j)=sinc(xx)+(sinc(xx-pi)+sinc(xx+pi))/2
         yout(j)=yout(j)*(1.0-(float(j)/(nhw+0.5))**2)**2 ! Apodize
         top=top+yout(j)
      end do
      top=2*top   ! Sum of the 2*NHW+1 operator values
c
c  Main loop:  Convolve interferogram with smoothing operator.
c  At the same time, find the peak AC interferogram value & location
      pinv=0.0
      pinl=1+nhw
      do i=1+nhw,nip-nhw
         smoo=yin(i)        ! Remember:  yout(0) = 1.0
         do j=1,nhw
            smoo=smoo+yout(j)*(yin(i-j)+yin(i+j))
         end do
         yout(i)=smoo/top
         yac=yin(i)-yout(i)
         if( abs(yac) .gt. abs(pinv) ) then
            pinv=yac
            pinl=i
         endif
      end do ! i=1+nhw,nip-nhw
c      write(*,*)'smooth_igram: pinl,pinv=',pinl,pinv

c  Extrapolate the smoothed interferogram into the end zones:
c  First NHW points (overwriting the no-longer-needed operator).
      smoo=yout(1+nhw)
      ds=yout(1+nhw)-yout(2+nhw)
      do i=nhw,1,-1
         smoo=smoo+ds
         yout(i)=smoo
      end do
c
c  Last NHW points.
      smoo=yout(nip-nhw)
      ds=yout(nip-nhw)-yout(nip-nhw-1)
      do i=nip-nhw+1,nip
         smoo=smoo+ds
         yout(i)=smoo
      end do

c  Fit a straight line across ZPD region of smoothed igram to
c  remove any artifact (e.g. dip) due to detector non-linearity.
c  The fit is weighted by the inverse of the AC interferogram
c  so that the points around ZPD have little impact.
       t0=0.0d0
       t1=0.0d0
       t2=0.0d0
       ty0=0.0d0
       ty1=0.0d0
       thresh=0.00004*pinv
       nwid=5*nhw ! half-width of ZPD region to be replaced by line
c  Max values of nwid= 5x600=3000
c
c  Don't let ZPD be too close to ends of interferograms
c  Such cases tend to be garbage scans anyway and can cause
c  very small OPD values in the runlog which cause GFIT problems
       if (pinl.le.nwid) then
           write(*,*) 'Warning: PINL < NWID'
           write(*,*) 'ZPD is too close to end of igram'
           write(*,*) 'Setting PINL to NWID+1'
           pinl=nwid+1
       endif
       if (pinl.gt.nip-nwid) then
           write(*,*) 'Warning: PINL > NIP-NWID'
           write(*,*) 'ZPD is too close to end of igram'
           write(*,*) 'Setting PINL to NIP-NWID'
           pinl=nip-nwid
       endif
c       if(pinl-nwid.lt.1)  pinl=1+nwid
c       if(pinl+nwid.gt.nip) pinl=nip-nwid
c
       do ii=-nwid,nwid
          yac=yin(pinl+ii)-yout(pinl+ii)
          wt=1/(1+(yac/thresh)**2)
          t0=t0+wt
          t1=t1+wt*ii
          t2=t2+wt*ii**2
          ty0=ty0+wt*yout(pinl+ii)
          ty1=ty1+wt*yout(pinl+ii)*ii
       end do
       det=t2*t0-t1**2  ! Determinant of 2x2 matix
       if(det.eq.0.0d0) then 
          write(*,*)' dc_correction: SLF matrix is rank-deficient'
          dclevel=ty0/t0
          grad=ty1/t2
       else
          dclevel=(t2*ty0-t1*ty1)/det
          grad=(t0*ty1-t1*ty0)/det
       endif

c  Compute fractional size of ZPD artifact. Negative values imply dip
       frzpda=dble(yout(pinl))/dclevel-1.

c  Modify smoothed igram to remove artifact/dip at ZPD. 
c  Use a linear combination of the fitted straight line
c  and the original smoothed interferogram: At ZPD use 100%
c  the fitted straight line; at the edges of the ZPD window,
c  use 100% smoothed igram.  This "windowing" avoids
c  discontinuities at pinl-nwid and pinl+nwid
      do ii=-nwid,nwid
         fwt=(1-(float(ii)/nwid)**2)**2
         yout(pinl+ii)=(1-fwt)*yout(pinl+ii)+fwt*(dclevel+ii*grad)
      end do
      return
      end

      real*8 function sinc(t)
      real*8 t,tt
      if (dabs(t).lt.0.2d0) then
         tt=t*t
         sinc=1.0d0-tt*(1.0d0-tt/20)/6  !  Approx good to 10^-8 at t=0.2
      else
         sinc=dsin(t)/t
      endif
      return
      end

