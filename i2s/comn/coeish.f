      subroutine compute_odd_even_igram_shift(lgev,ac_igram,
     & shbar,shrms)
c
c  Determines shift between even and odd interferogram subsets.
c  This can be non-zero due to uneven spacing (mis-sampling)
c  of the laser fringe crossings.
c
c  Inputs:
c     LGEV          I*4  Length of Ghost Estimation Vector
c   AC_IGRAM(LGEV)  R*4  Igram portion to be used for ghost estimation.
c
c  Output:
c    SHBAR          R*4  Odd/Even shift (as fraction of sample interval)
c    SHRMS          R*4  Standard Deviation of shift (same units)
c
c  Theoretical Basis:
c
c  FTIR spectrometers commonly use both the rising and the
c  falling laser mean crossings to trigger the sampling of
c  the IR interferograms. This doubles the free spectral range,
c  extending it up to the laser frequency, as compared with
c  sampling only on the rising or falling mean crossings.
c  If the laser mean level is biased, this results in an uneven
c  (short-long-short-long) sampling of the interferogram. All the
c  odd interferogram samples are separated by exactly the laser 
c  wavelength, and all the even samples too. But successive samples
c  are alternately shorter/longer than the nominal half-laser-wavelength.
c  The resulting periodic sampling error causes ghosts seperated
c  by the laser frequency which are then aliased into the frequency
c  range of interest, causing mixing of signals across the laser
c  half-frequency. 
c
c  Mis-sampling is estimated from the resulting rotation of the
c  phase curve. The non-ZPD-parity igram subset is rotated by
c  SHBAR+Pi/2 as compared with the ZPD-parity igram, the Pi/2
c  coming from the (forced) choice of an adjacent ZPD point.
c  And the phase curve computed from all the igram points is
c  rotated by SHBAR/2 with respect to that computed from the
c  ZPD-parity igram subset. This provides two ways to estimate
c  the phase rotation, both of which are used and their results
c  averaged.
c
c  Each point in the phase spectrum therefore provides an
c  estimate of the phase rotation, and hence SHBAR. The
c  outputted value of the phase rotation is a weighted
c  average of that from the individual points, and its
c  uncertainty (SHERR) is the standard error. 
c  The LSE-related phase rotation is usually very small
c  compared with the natural undulations of the phase curve,
c  and so care must be taken not to let any aspect of the
c  phase curve itself affect the phase differences (e.g.
c  2.Pi discontinuities in the two phase curves).
c
c  Since aliasing can bring spurious
c  signal into the domain of interest, carrying the wrong phase,
c  points affected by aliasing must be heavily de-weighted in
c  the averaging process. Fortunately, we can tell from the
c  full igram how much aliasing there will be at each point
c  in the igram-subsets used in the phase different calculation.
c  These aliasing-affected points can then be appropriately
c  de-weighted in the averaging of the phase differences.
c 
c    Assume the raw interferogram is band-limited, but mis-sampled.
c  If we construct a new interferogram comprising of just the
c  odd or even interferogram points, the sampling is now perfect,
c  but the interferogram is, in general, not band limited. [If it was,
c  then the interferogram size and data rate was a factor two
c  larger than necessary.]
c
c  Implementation Notes:
c    Input vector (AC_IGRAM) is unmodified.
c
c    ZPD location is at AC_IGRAM(LGEV/2+1) which is an 
c    odd point as far as THIS subroutine is concerned..
c
c    - j is used to index interferogram points
c    - k is used to index spectral points
c    - i = SQRT(-1)
c
c  The phase difference is evaluated directly from complex
c  spectra (not by subtracting the absolute phases) to avoid
c  the  2*Pi  discontinuity problem.
c  Since  TAN(x2-x1)= [TAN(x2)-TAN(x1)]/[1+TAN(x1)*TAN(x2)]
c  If we have two complex vectors [a1,b1] and [a2,b2] such that
c     TAN(x2)=b2/a2  and TAN(x1)=b1/a1
c  Then  TAN(x2-x1)= [b2/a2-b1/a1]/[1+(b1/a1)*(b2/a2)]
c  Hence TAN(x2-x1)= [a1.b2-a2.b1]/[a1.a2+b1.b2]
c  Hence x2-x1 = ATAN2( a1.b2-a2.b1, a1.a2+b1.b2)


      implicit none
      integer*4 j,ji,k,lgev,lgev2,lgev4,
     & icall,kflip,kg,iflip,
     & kr,ki,kf

      real*4
     & shbar,shrms,
     & shabar,sharms,
     & shbbar,shbrms

      real*4
c     & power1_all,
c     & aliasing0_was,da,dp,
     & phase0,
c     & phase1,
     & phasei,phasek,yy,
     & lasf,pi,
     & x,apo,
     & preda,predb,
     & ac_igram(lgev), ! input igram subset (unmodified)
     & wsall0(lgev),   ! working vector for fully-sampled igram
c     & wsall1(lgev),   ! working vector for fully-sampled igram
     & wsoe(lgev/2,2),  ! working array for under-sampled igrams
     & dnu

      real*8 dpa,cpa,dpb,cpb,
     & power0_all,tlower,
     & power_oei,power_oek,
     & aliasing0,
     & ar0,ai0,ar1,ai1,
     & den,wt,
     & delpa,delpb,
     & frdelpa,frdelpb,
     & toti2, tots0, totsi, totsk,
c     & totsa,totsd,
     & tw,tdela,tda2,tdelb,tdb2

      data icall/0/

      icall=icall+1

      pi=3.14159265
      tlower=-1.d0  ! Avoid compiler warnings
      lgev2=lgev/2
      lgev4=lgev/4
      lasf=15798.03
      dnu=lasf/lgev2

c  Initialize workspace vectors containing:
c  WSALL()  Full (odd+even) igram points
c  WSOE(*,1) odd-only (ZPD-parity) igram points
c  WSOE(*,2) even-only (non-ZPD-parity) igram points
c     write(54,*) 2,2
c     write(54,*) ' j wsall'
c     write(55,*) 2,3
c     write(55,*) ' k  wsoe1  wsoe2 '
      ji=1
      iflip=1
c     toti2=0.0d0
      do j=1,lgev
         x=float(j-1-lgev2)/lgev2
         apo=(1-x**2)**2
         yy=apo*ac_igram(j)
         wsall0(j)=yy
c         wsall1(j)=yy
         wsoe(ji,iflip)=2*yy !  Double WSOE igrams (half the points)
c        if(iflip.eq.2) write(55,*)ji,wsoe(ji,1),wsoe(ji,2)
         ji=ji+iflip-1
         iflip=3-iflip
c        toti2=toti2+yy**2
c        write(54,*) j,wsall0(j)
      end do

c  Convert igrams into complex spectra.
c  Rotate by LGEV/2, such that ZPD is at WSALL(1), then FFT
      call vrot(wsall0,lgev2,lgev)   ! FFT-pack
      call ffak(wsall0,lgev)         ! Fourier Analysis
c      call vrot(wsall1,1+lgev2,lgev) ! FFT-pack
c      call ffak(wsall1,lgev)         ! Fourier Analysis

c  Look at full-bandwidth spectrum from fully-sampled (odd+even)
c  igram to see where (and how much) aliasing will occur when
c  the odd-only and even-only interferogram subsets are FFT'd.
c
c  This spectrum is used subsequently to weight the phase
c  difference spectrum to mitigate the effects of aliasing
c  and 1/f noise while using as much unaliased spectral signal
c  as possible.
c
c  Determine whether signal is in upper or lower half of alias.
c      write(56,*) 2,7
c      write(56,*) ' ic  k   wav   amp0   ph0/pi   amp1   ph1/pi '
c      write(56,56)icall,0,0*dnu,abs(wsall0(1)),.0,abs(wsall1(1)),.0
c56    format(2i4,7f12.5)
      tots0=0.5*(dble(wsall0(1))**2+dble(wsall0(2))**2)  ! DC & Nyquist terms
      do k=1,lgev2-1
         kr=2*k+1  ! pointer to real element
         ki=2*k+2  ! pointer to imag element
         power0_all=dble(wsall0(kr))**2+dble(wsall0(ki))**2
c         power1_all=wsall1(kr)**2+wsall1(ki)**2
         phase0=atan2(wsall0(ki),wsall0(kr))
c         phase1=atan2(wsall1(ki),wsall1(kr))
         if(k.eq.lgev4) tlower=tots0+power0_all/2
         tots0=tots0+power0_all
c         write(56,56)icall,k,k*dnu,sqrt(power0_all),phase0/pi,
c     &   sqrt(power1_all),phase1/pi
      end do
c      write(56,*)icall,k,k*dnu,abs(wsall0(2)),.0,abs(wsall1(2)),.0

      if(tlower.gt.0.5*tots0) then
         kflip=1   !  signal is mainly in lower half of alias
c      write(*,*) 100*tlower/tots0,'% of Signal in lower half of alias'
      else
         kflip=-1  !  signal is mainly in upper half of alias
c      write(*,*)100*(1-tlower/tots0),'% Signal in upper half of alias'
      endif

c  FFT odd-only (WSOE(*,1)) and even-only (WSOE(*,2)) igrams into spectra
      call vrot(wsoe(1,1),lgev4,lgev2)   ! FFT-pack (ZPD at WSOE (1,1))
      call ffak(wsoe(1,1),lgev2)         ! Fast Fourier Analysis
      call vrot(wsoe(1,2),lgev4,lgev2)   ! FFT-pack (ZPD at WSOE (1,2))
      call ffak(wsoe(1,2),lgev2)         ! Fast Fourier Analysis

c  Compute phase differences between odd/even spectra
c  and find their weighted average and its uncertainty.
c     write(57,*)' 2 16 '
c     write(57,*) ' k   wavn     amp0      amp1      ampo     ampe   
c    &ph0/pi    ph1/pi    pho/pi    phe/pi   aliasing     wt     delpa/p
c    &i  delpb/pi     fra     frb'
         tw=0.0d0
         tdela=0.0d0
         tda2=0.0d0
         tdelb=0.0d0
         tdb2=0.0d0
         totsi=0.5*(dble(wsoe(1,1))**2+dble(wsoe(2,1))**2) ! Parseval: DC & Nyquist
         totsk=0.5*(dble(wsoe(1,2))**2+dble(wsoe(2,2))**2) ! Parseval: DC & Nyquist
c  If spectrum is aliased (kflip=-1) then kf=k
c  If spectrum is unaliased (kflip=1) then kg=k
         kf=(kflip+1)*lgev4  ! kf = 0 or lgev/2
         kg=(1-kflip)*lgev4  ! kg = lgev/2 or 0
         do k=1,lgev4-1
            kf=kf-kflip
            kg=kg+kflip
            kr=2*k+1  ! pointer to real element
            ki=2*k+2  ! pointer to imag element
            power0_all=dble(wsall0(2*kg+1))**2+dble(wsall0(2*kg+2))**2
c            power1_all=wsall1(2*kg+1)**2+wsall1(2*kg+2)**2
            phase0=atan2(wsall0(2*k+2),wsall0(2*k+1))
c            phase1=atan2(wsall1(2*k+2),wsall1(2*k+1))
            aliasing0=dble(wsall0(2*kf+1))**2+dble(wsall0(2*kf+2))**2
c
            ar0=wsoe(2*k+1,1)
            ai0=wsoe(2*k+2,1)
            ar1=wsoe(2*k+1,2)
            ai1=wsoe(2*k+2,2)
            cpa=ar1*ai0-ar0*ai1
            dpa=ar0*ar1+ai0*ai1
            cpb=wsall0(2*kg+1)*ai0-kflip*ar0*wsall0(2*kg+2)
            dpb=ar0*wsall0(2*kg+1)+kflip*ai0*wsall0(2*kg+2)
            delpa=datan2(cpa,dpa)/pi   ! Phase diff (Odd-only minus Even-only)
            delpb=datan2(cpb,dpb)/pi   ! Phase diff (Odd-only minus Full)

c  Predict phase difference in absense of LSE
            preda=kflip*(1-float(kf)/lgev2)   ! Odd-only minus Even-only
c            predb=0.25*(kflip-1)             ! Odd-only minus Full
            predb=0.

c  Correct n.Pi phase discontinuities by using predictions
            delpa=delpa-float(nint(delpa-preda))
            delpb=delpb-float(nint(delpb-predb))

c  Express as fraction of sampling interval
            frdelpa=abs(lgev2*delpa/kg)-1
            frdelpb=kflip*lgev*(delpb-predb)/kg

c  Define weights to suppress 1/f noise and frequencies
c  damaged by aliasing.
            power_oei=dble(wsoe(kr,1))**2+dble(wsoe(ki,1))**2
            power_oek=dble(wsoe(kr,2))**2+dble(wsoe(ki,2))**2
            totsi=totsi+power_oei
            totsk=totsk+power_oek
            phasei=datan2(ai0,ar0)
            phasek=datan2(ai1,ar1)
            wt=power0_all**2/(power0_all+90000*aliasing0)
            tw=tw+wt
            tdela=tdela+wt*frdelpa
            tda2=tda2+wt*frdelpa**2
            tdelb=tdelb+wt*frdelpb
            tdb2=tdb2+wt*frdelpb**2
c            write(57,57)k,kg*dnu,sqrt(power0_all),sqrt(power1_all),
c     &      sqrt(power_oei),sqrt(power_oek),
c     &      phase0/pi,phase1/pi,phasei/pi,phasek/pi,
c     &      sqrt(aliasing0),wt,delpa,delpb,frdelpa,frdelpb
c57          format(i3,f8.1,18f10.6)
         end do   !   do k=1,lgev4-1
         kf=kf-kflip
         kg=kg+kflip
         power0_all=dble(wsall0(2*kg+1))**2+dble(wsall0(2*kg+2))**2
c         power1_all=wsall1(2*kg+1)**2+wsall1(2*kg+2)**2
         phase0=atan2(wsall0(2*kg+2),wsall0(2*kg+1))
c         phase1=atan2(wsall1(2*kg+2),wsall1(2*kg+1))
         aliasing0=dble(wsall0(2*kf+1))**2+dble(wsall0(2*kf+2))**2
         wt=power0_all**2/(power0_all+90000*aliasing0)
c         write(57,57)lgev4,kg*dnu,sqrt(power0_all),sqrt(power1_all),
c     &   abs(wsoe(2,1)),abs(wsoe(2,2)),phase0/pi,phase1/pi,0.0,0.0,
c     &   sqrt(aliasing0),wt,0.0,0.0,0.0,0.0
         shabar=tdela/tw
         sharms=sqrt( tda2*tw - tdela**2 )/tw
         shbbar=tdelb/tw
         shbrms=sqrt( tdb2*tw - tdelb**2 )/tw
         den= 1/sharms**2 + 1/shbrms**2
         shbar=(shabar/sharms**2+shbbar/shbrms**2)/den
         shrms=1/dsqrt(den)
c         write(*,'(6f11.7)') shabar,sharms,shbbar,shbrms,shbar,shrms

      return
      end
