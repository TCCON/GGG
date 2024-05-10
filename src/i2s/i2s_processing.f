      subroutine i2s_processing(mip,mlong,verbose,idir,pinl,ichan,lpco,
     & pco_thresh,nburst,phasepath,ac_igram,nip,nlong,nshort,
     & nffthr,nside,
     & izpd,
     & zpa,shbar_out,sherr_out,lsemode,fpilha)

c  Processes a raw interferogram into a spectrum. This
c  includes phase correction, ghost correction, and FFT.
c
c  Input:
c    mip           I*4    Maximum number of input points
c    mlong         I*4    Maximum half-FFT size
c    verbose       I*4    Level of verbosity for displayed messages
c    idir          I*4    Run direction 1=FWD ; (-1 or 0)=REV
c    pinl          I*4    Peak Interferogram Location (AC)
c    ichan         I*4    Channel (InGaAs=2; Si=1)
c    lpco          I*4    Length of phase correction operator
c    pco_thresh    R*8    Amplitude threshold for phase correction
c    nburst        I*4    Search width for center burst
c    phasepath     C*(*)  Where to write the phase curve
c    nip           I*4    Number of igram points
c    lsemode       I*4    The type of LSE correction to be done.
c
c  Input/output:
c    ac_igram(mip) R*4    Buffer containing Input_igram/Output_spectrum
c
c  Output:
c    nlong         I*4    # points on long side of phase-corrected igram
c    nshort        I*4    # points on short side of phase-corrected igram
c    nffthr        I*4    FFT size for igram-to-spectrum conversion
c    nside         I*4    Single-sided igram nside=1; double-sided nside=2
c    izpd          I*4    Point index (location) of ZPD
c    zpa           R*8    ZPD interferogram amplitude (phase-corrected)
c    shbar_out     R*4    Applied LSE shift
c    sherr_out     R*4    Applied LSE shift uncertainty
c    fpilha        R*4    Fraction of Power in Lower Half of Alias
c
      implicit none

      integer*4
     & mip,          ! Subroutine input argument (see above)
     & mlong,nlong,  ! Subroutine input argument (see above)
     & nshort,
     & nlimit,       ! = nshort (DS case); = nlong (SS case)
     & verbose,      ! Subroutine input argument (see above)
     & idir,         ! Subroutine input argument (see above)
     & nip,          ! Number of igram points
     & izpd,         ! Subroutine output argument (see above)
     & izpd_si,      ! Subroutine output argument (see above)
     & nffthr,       ! Half size of the high-res FFT
     & nfftpc,       ! Half-size of FFT needed for phase correction
     & nrip2bpc,     ! Number of raw igram points to be phase corrected
     & zff,          ! Zero-Fill-Factor
     & pinl,         ! From m4head.inc
     & nside,        ! 1 = single-sided; 2 = double-sided
     & ichan,        !
     & margin,       ! =nburst+lpco/2
     & conv_method,
     & ifurip,       ! Index of first used raw igram point
     & ilurip,       ! Index of last used raw igram point
     & nurip,        ! Number of used raw igram points = ilurip-ifurip+1
     & nhw,
     & lgev,
     & icall,
     & i,j,          ! indices
     & mpco,lpco,l2, ! Length of phase correction operator
     & lsemode,      ! Laser sampling error type (none,slave,master,Hase,other)
     & nburst 

      parameter (
     & mpco=1024*64, ! Maximum Length of Phase Correction Operator
     & lgev=1024,    ! Length of Ghost Estimation Vector
     & conv_method=2) ! 1=Slow, 2=fast convolution

      real*4 temp, zero,shpmo,
     & y_nyquist,
     & shbar,sherr,fpilha,
     & shbar_out,sherr_out,
     & ac_igram(mip),  ! Subroutine input/output argument (see above)
     & bufn(mip/2),
     & best,
     & pco(mpco)    ! Phase Correction Operator

      real*8
     & dfr,        ! ZPD at IZPD+DFR
     & zpa,        ! ZPD interferogram amplitude (phase-corrected)
     & pco_thresh, ! Amplitude threshold for phase correction
     & zpdl        ! From m4head.inc

      character
     & phasepath*(*) ! path to write phase file

      save shbar,sherr,izpd_si
      data icall/0/

      write(*,*) 'Entering i2s_processing'
      icall=icall+1
      zero=0.0
      if(lpco.gt.mpco) stop 'Error in i2s_processing: LPCO > MPCO'

c  REV interferogram are reversed in-place
      if(idir.le.0) then
         do i=1,nip/2
            temp=ac_igram(i)
            ac_igram(i)=ac_igram(nip+1-i)
            ac_igram(nip+1-i)=temp
         end do
         pinl=nip+1-pinl  ! ZPD moves too
      endif

c  Check that the ZPD point is more than NBURST+LPCO/2 from the
c  ends of the interferogram (to avoid an array-bound error later).
      margin=nburst+lpco/2
      l2=lpco/2
      if(pinl.lt.1+margin) then
         write(*,*) 'pinl, margin = ',pinl,margin
         stop 'PINL < 1+Nburst+LPCO/2'
      endif
      if(pinl.gt.nip-margin) then
         write(*,*) 'pinl, nip-margin = ',pinl,nip-margin
         stop 'PINL > NIP-Nburst-LPCO/2'
      endif

c  Best ZPD point is not necessarily the one having the peak AC
c  interferogram amplitude, expecially if igram is chirped.
c  Find interferogram point closest to the maximum in symmetry.
      call bestzpd(ac_igram(pinl-margin),lpco,nburst,best)
      zpdl=dble(pinl)+dble(best)
      izpd=nint(zpdl)
      dfr=zpdl-izpd
      dfr=0.0d0
c      izpd=pinl   !  Old way of defining ZPD
      if(verbose.ge.3) write(*,*)'i2s_processing:'//
     & ' Best ZPD point: pinl,izpd,zpdl = ',pinl,izpd,zpdl

c  Ghost Correction
      if(ichan.eq.2) then
         call compute_odd_even_igram_shift(lgev,ac_igram(izpd-lgev/2),
     &   shbar,sherr,fpilha)
         izpd_si=izpd
      endif

c  Only re-sample igram if LSEMODE=2 and nearly all spectral power confined
c  to one half of alias. Shift & resample non-ZPD-parity interferogram points
c  by shbar/2.  Points with ZPD-parity (odd/evenness) as IZPD are unchanged.
c  SHBAR is halved because the points in BUFN have a spacing of 2 and
c  the SHARS subroutine treats SHBAR as the fraction of the BUFN spacing.
      if(lsemode.eq.2 ) then ! Re-sample non-ZPD-parity igram points
         if(fpilha*(1.-fpilha).gt.0.001 ) write(*,*) 'Warning: '//
     &   'FPILHA  > 0.001  SHBAR estimate is likely poor'
         call vmov(ac_igram,1,bufn,1,nip)        ! Copy igram points to BUFN
         shpmo=-shbar*(-1)**izpd_si
         nhw=60
         if(ichan.eq.2) then  ! Si
            call resample_ifg(nip,bufn,ac_igram,1,shpmo,nhw)
         else                ! InGaAs
            call shars(1,nip,nhw,shpmo,bufn,ac_igram)
         endif
         shbar_out=shbar
         sherr_out=sherr
c         write(*,*) ' LSE correction performed: ',lsemode,fpilha
      else
         nhw=0
         shbar_out=0.0
         sherr_out=0.0
c         write(*,*) ' LSE correction skipped: ',lsemode,fpilha
      endif

cc  Recompute shifts to check that they have actually gone to zero.
c      if(ichan.eq.2) then
c        call compute_odd_even_igram_shift(lgev,ac_igram(izpd-lgev/2),
c     &  shbar2,sherr2,fpilha2)
cc        write(*,*) 'ichan,shbar2=',ichan,shbar2/2,sherr2/2,idir,izpd
c      endif

c  Phase Correction
c
c  Copy igram center-burst into PCO in FFT-packed format
      call vmov(ac_igram(izpd),1,pco,1,l2)
      call vmov(ac_igram(izpd-l2),1,pco(l2+1),1,l2)

c  Calculate phase correction operator
      call phase_corr_oper(lpco,pco,dfr,pco_thresh,phasepath)

c  Set FFT size to accept the full interferogram without exceeding
c  the limit 'mlong' obtained from the input file.
c      nlong=min(nip-izpd-l2,mlong)
c      nshort=min(izpd-1,nip-izpd)
      nlong=max(izpd-1,nip-izpd)-l2
      nshort=min(izpd-1,nip-izpd)-l2
      if(nlong.gt.mlong) then
         write(*,*) 'Warning: Truncating igram from NLONG = ',nlong,
     & ' to MLONG = ',mlong
         nlong=mlong
      endif
      if(nshort.gt.mlong) then
         write(*,*) 'Warning: Truncating igram from NSHORT = ',nshort,
     & ' to MLONG = ',mlong
         nshort=mlong
      endif
c     write(*,*) 'i2s_processing: nip,izpd,l2,mlong=',nip,izpd,l2,mlong

      if(iabs(izpd-nip/2) .gt. nip/6) then
         nside=1  ! Single-sided FWD
         ifurip=nshort+2  ! Index of first used raw igram point
         nurip=nlong+2*l2
c         write(*,*) ' Single-sided',nside
         nlimit=nlong
      else
         nside=2  ! Double-sided FWD
         ifurip=1  ! index of first used raw igram point
         nurip=2*(nshort+l2)+1
c         write(*,*) ' Double-sided',nside
         nlimit=nshort
      endif
      ilurip=ifurip+nurip-1  ! Index of last used raw igram poiont

c  Determine FFT-size for phase correction
      nfftpc=1
      do while (nfftpc.lt.nurip)
         nfftpc=nfftpc*2
      end do
      if(nfftpc.gt.mip) then
         write(*,*)'nfftpc, mip=',nfftpc,mip
         stop 'Increase parameter MIP'
      endif
c     write(*,*)'i2s_processing: mlong,nlong,nrip2bpc,nfftpc = ',
c    & mlong,nlong,nrip2bpc,nfftpc

c  2018-12-22
c  Phase Correction will use the raw igram points from:
c  Izpd-L2+1 to NIP = Nlong+Izpd, a total of Nlong+L2 points.
c  Points outside this range are set to zero if conv_method=2
c  as a precaution; in case they contain a very large value
c  that will get sucked into the FFT and spread its rounding
c  error about.
c  Since the convolution origin is shifted from 1 to IZPD-L2+1, to
c  minimize the required size of NFFTPC, the igram must be zeroed from
c  Nlong+IZPD+1 to NFFTPC+IZPD-L2, a total of NFFTPC-L2-Nlong points

c  Phase Correction will utilize igram points: Izpd-L2 to Nlong+L2+Izpd
c  No need to zero-fill from 1 to Izpd-L2-1, it won't be seen by CONVOL
c  Zero-fill the uninitialized section of ac_igram beyond Nlong+L2+Izpd
c  Since the convolution origin is shifted from 1 to IZPD-L2, to
c  minimize the required size of NFFTPC, the igram must be zeroed from
c  Nlong+IZPD+L2+1 to NFFTPC+IZPD-L2, a total of NFFTPC-LPCO-Nlong points

c 2018-03-13
c  Phase Correction uses igram points: Izpd-L2+1 to NIP-NHW = Nlong+Izpd-NHW
c  No need to zero-fill from 1 to Izpd-L2, it won't be seen by CONVOL.
c  Zero-fill the uninitialized section of ac_igram beyond NIP=Nlong+Izpd-NHW
c  Igram zeroed from  Nlong+Izpd-NHW+1 to NFFTPC+Izpd-L2, a total of
c  NFFTPC-L2-Nlong+NHW  points

c      write(*,*)'zeroing: nip+1  to  nfftpc =',nxxlong+izpd+l2+1,
c     & izpd+nfftpc
c  Zero the parts of the igram that won't be needed as a pre-caution
      if(conv_method.eq.2) then
       call vmov(zero,0,ac_igram(1),1,ifurip-1)  ! before first used point
c      call vmov(zero,0,ac_igram(ilurip+1),1,nfftpc-ilurip) ! after last
       call vmov(zero,0,ac_igram(ilurip+1),1,nfftpc+ifurip-ilurip-1) ! after last
      endif

c  Perform in-place convolution of raw igram in AC_IGRAM with PCO
c    After convolution, valid igram extends from IZPD to izpd+nlong
c      call convol(pco,lpco,ac_igram(1),nfftpc)
      call convolve(conv_method,nfftpc,ac_igram(ifurip),bufn,
     & lpco,pco)
      zpa=dble(ac_igram(izpd))  !  ZPD igram amplitude (phase-corrected)

c      write(*,*)'izpd, zpa =',izpd,zpa 
c      write(78,*)' 2 2 '
c      write(78,*)' i y '
c      do i=1,1024*2048
c        write(78,*)i,ac_igram(i)
c      end do
c      close(78)

c  High-Res FFT size for FFAK
      zff=1  ! Zero fill factor
      nffthr=1
      do while (nffthr.lt.nside*nlimit)
         nffthr=nffthr*2
      enddo
      if(nffthr.gt.mip/nside) then
         write(*,*)'nffthr, mip=',nffthr,mip
         stop 'i2s_processing: Increase parameter MIP'
      endif
      nffthr=nffthr*zff
c     write(*,*)'i2s_processing:  NFFThr,nside = ',nffthr,nside

      if(nside.ge.2) then ! double-sided
c         write(*,*) 'double-sided igram'
         call vrot(ac_igram,izpd-1,nffthr) ! circular shift brings ZPD to point 1
         call vmov(zero,0,ac_igram(nlong+1),1,nffthr-nlong-izpd+l2) !  Zero-fill
      else               ! single-sided
c         write(*,*) 'single-sided igram'
         call vmov(ac_igram(izpd),1,ac_igram(1),1,nlong) ! bring ZPD from Izpd to 1
c  Zero-fill unused portions of AC_IGRAM beyond long side of interferogram,
c  including the NFFT+1'st point which otherwise would be ioverlooked.
         call vmov(zero,0,ac_igram(nlong+1),1,nffthr-nlong+1) ! Zero-fill

c  Unfold long side of interferogram, together with appended DC-filled points
         j=2*nffthr
         do i=2,nffthr
            ac_igram(j)=ac_igram(i)
            j=j-1
         end do
      endif  !  if(nside.ge.2) then ! double-sided

c      write(66,*) 2,3
c      write(66,*) ' i  y_left y_right'
c      write(66,*) 1, ac_igram(1), ac_igram(1)
c      do i=2,nffthr/2
c         write(66,*) i, ac_igram(i), ac_igram(2*nffthr+2-i)
c      end do

c  Real-to-complex FFT of double-sided, symmetrical interferogram
      call ffak(ac_igram,2*nffthr/nside)

c      write(68,*) 2,2
c      write(68,*) ' i  real '
c      write(68,*) 0, 0.125*ac_igram(1) ! DC term
c      do i=1,nffthr/nside-1
c         write(68,*) i, 0.125*ac_igram(2*i+1)
c      end do
c      write(68,*) i, 0.125*ac_igram(2) ! Nyquist

c  Perform compaction, i.e. overwrite imaginary values, which are
c  all zero anyway in SS case, with the real values. This halves the
c  length of the vector. Also normalize spectral values to match OPUS.
c  Need to pull out the Nyquist term first (since it is located in the
c  first imaginary location) in case the saved spectrum will extend
c  to Nyquist.
      y_nyquist=ac_igram(2) !  GCT 20190310
      call vmul(ac_igram,2,0.125,0,ac_igram,1,nffthr/nside)
      ac_igram(nffthr/nside+1)=0.125*y_nyquist  !  GCT 20190310
c  The statement above, that scales the Nyquist amplitude, although
c  correct, inexplicably changes the values of most of the spectrum.
c  This problem appears to have been corrected by the vmov indexing
c  change on line 261 in July 2020.


c      write(*,*) 'Exiting i2s_processing'
      return
      end

c================================================
      subroutine bestzpd(ac_igram,np,nburst,best)
c  Scans the interferogram from -NBURST to NBURST and finds the point or half-
c  point having the best symmetry. Then fits a parabola to the symmetry at
c  this point (or half point) and its two immediate neighbors to determine the
c  fractional interferogram location at which the parabola was a maximum.
c  This location is then returned as the best symmetry
c
c  INPUTS:
c     AC_IGRAM      the raw interferogram
c     NP           number of points over which to compute symmetry
c     NBURST       number of interferogram points to be tested
c
c  OUTPUTS:
c     BEST         fractional interferogram point having highest symmetry
c
c  Note that the symmetry is evaluated every half point, and then the 3
c  half points bracketing the point of maximum symmetry are fitted by a
c  parabola in order to determine the fractional location of peak symmetry.
c  The best ZPD point to choose would usually be NINT(BEST).

      implicit none

      integer*4 i,np,nburst
      real*4 ac_igram(-nburst-np/2:nburst+np/2),
     & symmi,symmp,symiw,sympw,best,smax,eps,denom
c
      eps=1.e-37
      smax=-999.0
      best=0.0
      symiw=0.0
      sympw=0.0
      do i=-nburst,+nburst
         call symmetry(ac_igram(i-np/2),symmi,symmp,np)
         if(sympw.gt.smax) then
            smax=sympw
            denom=eps+4.0*abs(2.0*sympw-symiw-symmi)
            best=float(i)-0.5+(-symiw+symmi)/denom
         endif
         if(symmi.gt.smax) then
            smax=symmi
            denom=eps+4.0*abs(2.0*symmi-sympw-symmp)
            best=float(i)+(-sympw+symmp)/denom
         endif
         symiw=symmi
         sympw=symmp
      enddo
      return
      end
