      function compute_snr(nmp,y,nsmoo,fors)
c  Estimates the Signal-to-Noise Ratio of an absorption spectrum.
c  Signal is estimated from the peak smoothed value.  Noise is
c  estimated from the rms differences between sucessive smoothed
c  points.
c
c  Inputs: 
c     NMP    I*4  Number of points in input vector.
c    Y(NMP)  R*4  Input vector containing spectrum (unchanged on output)
c     NSMOO  I*4  Half-Width of trapezoidal smoothing operator.
c     FORS   R*4  Fraction Of RMS Signal above which spectrum is de-weighted.
c
c  Outputs:
c      SNR  R*4  Signal-to-Noise Ratio
c
c
c  Conceptually, the spectrum is smoothed. We look at the rms diffs
c  between the high-res and smoothed spectra in regions where the
c  signal is close to zero.  If the smoothed signal is not close to 
c  zero, then there are absorption lines possibly present which will
c  drive up the RMS. The problem with this approach are noise spikes
c  at low wavenumbers where the smoothed spectrum is still close to
c  zero but the differences are huge. Ditto for 1/f noise, but this
c  can be easily suppressed by applying a f-dependent weighting to
c  the spectrum. Perhaps an abs(y-ysmoo) weighting would diminish
c  the effect of the noise spikes.
c
c  Investigate  Sum |y-ysmoo|/N  versus Sqrt (Y-Ysmoo)^2 /
c
c
c  SIGNAL.
c  The signal is computed from the smoothed spectrum. The smoothing
c  prevents noise spikes (especially at low wavenumber) from being
c  mistaken for signal.  The signal level is defined as
c     abs(ymax+ymin)
c  In a good spectrum, ymin should be slightly negative, whereas
c  ymax is strongly positive.  In a bad spectrum with ymax and ymin
c  equal and opposite (e.g. phase error) the signal will be ~0.
c
c  NOISE.
c  The noise estimation is more complicated. The noise is estimated
c  from the differences between successive smoothed points, rather
c  than their deviation from zero.  This has the advantage 
c  of being less sensitive to zero level offsets. Subroutine assumes
c  that if signal is close to zero, this is a good place to estimate
c  the noise level, since such regions are likely blacked out by gas
c  absorptions or be out-of-band. So the near-zero-signal points are
c  highly weighted in computing the noise. Conversely, if the signal
c  is above average, this is a bad place to estimate the noise level
c  due to the possible presence of sharp absorption lines. So these
c  points are de-weighted.
c
c  Another type of weighting is done depending on the index of 
c  the points within the input vector (Y). Points at the ends are
c  de-weighted in comparison with those in the middle, to counter
c  the effects of 1/f noise, low-frequency spikes, and aliasing.
c
c SMOOTHING.
c  Spectrum is smoothed to prevent ringing from the sinc from
c  being interpreted as noise, or to prevent noise spikes from
c  being interpreted as signal. For smoothing we use a trapezoid
c  filter whose end points are half the weight of the inner points.
c  If the half-width of the filter is N, then the weights are
c      w(-N)   = 1/4N
c      w(-N+1) = 1/2N
c      w(-N+2) = 1/2N
c        .        .
c        .        .
c      w(0)    = 1/2N
c        .        .
c        .        .
c      w(N-2)  = 1/2N
c      w(N-1)  = 1/2N
c      w(N)    = 1/4N
c  The total number of operator points is 2N+1
c  The sum of all the weights is (2N+1)*(1/2N)-2*(1/4N) = 1
c
c     Ysmoo(i) = Sum{w(k).y(i+k)} k=-N,+N
c     Ysmoo(i) = y(i-N)/4N + y(i+N)/4N + Sum{y(i+k)/2N} k=-N+1,+N-1
c     Ysmoo(i+1) = y(i+1-N)/4N + y(i+1+N)/4N + Sum{y(i+1+k)/2N} k=-N+1,+N-1
c
c     Ysmoo(i) = y(i-N)/4N + y(i+N)/4N + y(i-N+1)/2N + Sum{y(i+k)/2N} k=-N+2,+N-1
c     Ysmoo(i+1) = y(i+1-N)/4N + y(i+1+N)/4N + Sum{y(i+k)/2N} k=-N+2,+N
c     Ysmoo(i+1) = y(i+1-N)/4N + y(i+1+N)/4N + y(i+N)/2N + Sum{y(i+k)/2N} k=-N+2,+N-1
c
c  The Sum{} terms cancel when subtracting Ysmoo(i) from Ysmoo(i+1) yielding
c     Ysmoo(i+1)-Ysmoo(i) = [ - y(i-N) - y(i-N+1) + y(i+1+N) + y(i+N) ]/4N
c     4N.Ysmoo(i+1)-4N.Ysmoo(i)) = - y(i-N) - y(i-N+1) + y(i+1+N) + y(i+N)
c
c  So 4N.Ysmoo(i+1) can be computed from 4N.Ysmoo(i) with just four
c  operations: two additions and two subtractions. This is independent
c  of the operator length, N, since the contributions of the inner
c  points of the operator are unchanged.
c
c  4N.Diff = 4N.Ysmoo(i+1) - 4N.Ysmoo(i)
c          = -y(i-n) - y(i-n+1) + y(i+n) + y(i+n+1) 
c
c  4N.Ybar = (4N.Ysmoo(i+1)+4N.Ysmoo(i))/2 
c          =  4N.Ysmoo(i) + 4N.Diff/2
c
c  The trapezoidal smoothing reduces rms deviation by a factor
c  Sqrt((4N-1)/8N^2). For the N=1 case of a triangular operator
c  with 3 points of weights [0.25, 0.50, 0.25] the reduction is
c  Sqrt(3/8)=0.6124. For large N this expression tends to
c  1/Sqrt(2N) which is the square root of the operator FWHM.
c  If the weights were all equal (i.e. rectangular operator),
c  then the noise would be reduced by a factor 1/Sqrt(2N+1),
c  which for the N=1 case would be  1/Sqrt(3) = 0.577.
c  
c  SNR is calculated as the peak smoothed signal divided by the
c  noise level. For random noise the point-to-point variation is
c  1.414 times larger than the rms noise level. But after smoothing
c  the noise is no longer random, and so this 1.414 factor no
c  longer applies.
c
c  The table below shows the noise properties of random noise that
c  has been evaluated in two ways: 1) RMS deviation from the mean
c  or 2) RMS deviation from the previous point. This has been done
c  for various HWHM widths (NHW) of the trapezoid smoothing filter.
c
c   NHW       RMS            RMS          Ratio
c            Mean           Diff           M/D
c----------------------------------------------
c    0       1.000          1.414         0.707
c    1       0.612          0.500         1.225
c    2       0.468          0.250         1.871
c    3       0.391          0.167         2.347
c    4       0.342          0.125         2.739
c    5       0.308          0.100         3.082 
c    .         .              .             .
c    .         .              .             .
c    N   Sqrt(2N-0.5)/2N    1/2N      Sqrt(2N-0.5)
c-----------------------------------------------
c
c  For N > 0, the RMS deviation from the mean decreases with N
c  according to the expression
c      Sqrt(2N-0.5)/(2N)
c  For N > 0 the RMS deviation between sucessive points decreases
c  with N  according to the expression
c      1/(2N) 
c  The ratio of these expressions is
c     Sqrt(2N-0.5)
c
c  So to convert RMS difference between successive smoothed points
c  to RMS deviation from the mean in the raw unsmoothed data,
c  multiply by 2N. This expression is specific to the trapezoidal
c  smoothing operator with the end points having half the value of
c  the inner points.
c
c  The computed SNR value is scale-invariant: if the spectral values
c  were multiplied by a constant, even a negative constant, the
c  computed SNR value wouldn't change. The computed SNR value is also
c  flip-invariant: if the spectrum is reversed, the SNR is unchanged.
c
      implicit none
      integer*4 i,
     & imin,imax,
     & icall,kcall,
     & nmp,nsmoo
      real*4 y(nmp),ymin,ymax,compute_snr,enoise,diff4n,ymean,yrms,
     & fors,a2,xi,xbar,xwid,dnu,pi
      real*8 tot0,tot1,tot2,tot,wt,tw,twd1,twd2,
     & ysnex4n,ysmoo4n,yswas4n,sf
      parameter (pi=3.14159265)
      logical debug

      data icall/1/

      kcall=258
      debug=.false.
      debug=.true.
      sf=(0.5/nmp)**2
      dnu=0.007533
      dnu=1.0/nmp

c  Compute first two moments of spectral distribution (mean & width)
      tot0=0.0d0
      tot1=0.0d0
      tot2=0.0d0
      do i=1,nmp
         xi=float(i)
         tot0=tot0+y(i)
         tot1=tot1+y(i)*xi
         tot2=tot2+y(i)*xi**2
      end do
      xbar=tot1/tot0
      xwid=sqrt(tot2/tot0-xbar**2)
c     write(*,*)'xbar,xwid=',dnu*xbar,dnu*xwid
c     write(*,*)'ybar,yib=',tot0/nmp,tot0/xwid/sqrt(2*pi)

c  Ignoring DC term (i=1), compute Mean & RMS spectral values.
      tot1=0.5d0*y(nsmoo)            ! first point is given half-value
      tot2=0.5d0*y(nsmoo)**2         ! first point is given half-value
      do i=nsmoo+1,nmp-nsmoo-1
         tot1=tot1+y(i)          ! inner points are given full value
         tot2=tot2+y(i)**2       ! inner points are given full value
      end do
      tot1=tot1+0.5*y(nmp-nsmoo)     ! last point is given half-value
      tot2=tot2+0.5*y(nmp-nsmoo)**2  ! last point is given half-value
      ymean=tot1/(nmp-2*nsmoo)
      yrms=sqrt(tot2/(nmp-2*nsmoo))
      if(debug)write(*,*)'ymean,yrms = ',ymean,yrms
      a2=(4*nsmoo*fors*yrms)**2

c      if(icall.eq.kcall) write(67,*) 2,4
c      if(icall.eq.kcall) write(67,*) 'i y ysmoo wt'

c  Compute first smoothed point
      tot=0.5d0*y(1)   ! first point is half-value
      do i=2,2*nsmoo
         tot=tot+y(i)  ! interior points are full value
      end do
      tot=tot+0.5*y(2*nsmoo+1) ! last point is half-value
      yswas4n=2*tot    ! Normalize to 4N area
      i=nsmoo+1
      wt=sf*(i-1)*float(nmp-i)/(a2+yswas4n**2)
c      if(icall.eq.kcall) write(67,*) i,y(i),yswas4n/4/nsmoo,wt
      imin=i
      imax=i
      ymin=yswas4n
      ymax=yswas4n

c  Compute second smoothed point
      i=nsmoo+2
      ysmoo4n=yswas4n-y(1)-y(2)+y(2*nsmoo+1)+y(2*nsmoo+2)
      wt=sf*(i-1)*float(nmp-i)/(a2+yswas4n**2)
c      write(67,*) i,y(i),ysmoo4n/4/nsmoo,wt
      if(ysmoo4n.lt.ymin) then
         ymin=ysmoo4n
         imin=i
      endif
      if(ysmoo4n.gt.ymax) then
         ymax=ysmoo4n
         imax=i
      endif

      tw=0.0d0
      twd1=0.0d0
      twd2=0.0d0
      do i=nsmoo+3,nmp-nsmoo-1
         diff4n = -y(i-nsmoo-1)-y(i-nsmoo)+y(i+nsmoo-1)+y(i+nsmoo)
         ysnex4n = ysmoo4n + diff4n
         wt=sf*(i-1)*float(nmp-i)/
     &  (a2+yswas4n**2+ysmoo4n**2+ysnex4n**2)
c         if(icall.eq.kcall) write(67,*) i-1,y(i-1),ysmoo4n/4/nsmoo,wt
         twd1=twd1 + wt*abs(diff4n)
         twd2=twd2 + wt*diff4n**2
         tw = tw + wt

         if(ysmoo4n.lt.ymin) then
            ymin=ysmoo4n
            imin=i-1
         endif
         if(ysmoo4n.gt.ymax) then
            ymax=ysmoo4n
            imax=i-1
         endif

         yswas4n=ysmoo4n
         ysmoo4n=ysnex4n
      end do
c         if(icall.eq.kcall)write(67,*) i-1,y(i-1),ysmoo4n/4/nsmoo,wt
      if(debug) write(*,*)'compute_snr: imin,imax,n=',imin,imax,nmp
      if(debug) write(*,*)'compute_snr: ymin,ymean,ymax=',ymin/4/nsmoo,
     & ymean,ymax/4/nsmoo
c      enoise=2*nsmoo*sqrt(twd2/tw)
c  1.25 converts MAD to RMS deviation
      enoise=2*nsmoo*1.25*twd1/tw
      compute_snr=abs(ymax+ymin)/enoise
      if(debug)write(*,*)'comp_snr: sqrt(tdw/tw),enoise=',sqrt(twd2/tw),
     & enoise
c      close(67)
c     write(*,*) 'compute_snr: icall = ',icall
      icall=icall+1
      return
      end
