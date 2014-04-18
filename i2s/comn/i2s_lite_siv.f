      subroutine i2s_lite_siv (mip,mlfft,verbose,idir,pinl,ichan,lpco,
     & pco_thresh,phasepath,ac_igram,nip,izpd,zpa,
     & shbar_out,sherr_out,lsemode)
c
c  Input:
c    mip           I*4    Maximum number of input points
c    mlfft         I*4    Maximum log base 2 of half-FFT size
c    verbose       I*4    Level of verbosity for displayed messages
c    idir          I*4    Run direction 1=FWD ; (-1 or 0)=REV
c    pinl          I*4    Peak Interferogram Location (AC)
c    lpco          I*4    Length of phase correction operator
c    pco_thresh    R*8    Intensity threshold for phase correction
c
c  Input/output:
c    ac_igram(mip) R*4  Input interferogram/Output spectrum
c    nip           I*4    Number of input and output points
c
c  Output:
c    izpd          I*4    Point index (location) of ZPD
c    zpa           R*8    ZPD interferogram amplitude (phase-corrected)
c
      implicit none

      integer*4
     & mip,        ! Subroutine input argument (see above)
     & mlfft,      ! Subroutine input argument (see above)
     & verbose,    ! Subroutine input argument (see above)
     & idir,       ! Subroutine input argument (see above)
     & nip,        ! Subroutine input/output argument (see above)
     & izpd,       ! Subroutine output argument (see above)
     & izpd_si,    ! Subroutine output argument (see above)
     & nfft,       ! Half size of the high-res FFT
     & nfft2,      ! 
     & pinl,       ! From m4head.inc
     & fnbc,       ! Integer function First Non-Blank Character in string
     & indexa,     ! General loop index
     & ichan,      !
c    & koe,
     & kip,
     & lgev,
     & icall,nzero,
cgct     & maxp,
     & i,j,        ! indices
     & mpco,lpco,  ! Length of phase correction operator
     & lsemode,    ! Laser sampling error type (none,slave,master,Hase,other)
     & nburst 

      parameter (
     & mpco=1024*64,! Maximum Length of Phase Correction Operator
     & lgev=1024,   ! Length of Ghost Estimation Vector
     & nburst=15)  !

      real*4 temp, zero,
     & shbar,sherr,
     & shbar_out,sherr_out,
     & shbar2,sherr2,
     & ac_igram(mip),  ! Subroutine input/output argument (see above)
     & bufn(mip/2),
     & bufm(mip/2),
     & best,
     & cr(mpco)    ! Phase Correction Operator

      real*8
     & zpa,        ! ZPD interferogram amplitude (phase-corrected)
     & pco_thresh, ! phase correction intensity threshold
     & zpdl        ! From m4head.inc

      character
     & phasepath*(*), ! path to write phase file
     & stringa*11  ! String used to format integer display

      save shbar,sherr,izpd_si
      data icall/0/

      icall=icall+1
      zero=0.0
      if(lpco.gt.mpco) stop 'Error in i2s-lite-siv: LPCO > MPCO'

c  Set FFT size to accept the full interferogram without exceeding
c  the log-base-2 limit 'mlfft' obtained from the input file.
      nfft=1
      indexa=0
      do while((nfft.lt.nip).and.(indexa.lt.mlfft))
        nfft=nfft*2
        indexa=indexa+1
      enddo
      if(verbose.ge.4) write(*,*)'NIP, NFFT =',nip,nfft
      if((nfft*2).gt.mip) then
        write(stringa,'(i11)') (nfft*2)
        write(*,'(2a)')'Error: increase parameter MIP to ',
     &   stringa(fnbc(stringa):11)
        stop
      endif

c  Reverse interferogram are reversed in-place
      if(idir.le.0) then
ct        rund='R '
c        call vswap(ac_igram,1,ac_igram(nip),-1,nip/2)
        do i=1,nip/2
           temp=ac_igram(i)
           ac_igram(i)=ac_igram(nip+1-i)
           ac_igram(nip+1-i)=temp
        end do
        pinl=1+nip-pinl
      else
ct        rund='F '
      endif

c  KIP is the number of interferogram points that will be used
      kip=nip

cgctc  Limit size of interferogram to maxp points
cgct      maxp=1420200
cgct      if(kip.gt.pinl+maxp) then
cgct         write(*,*) 'Truncating interferogram to', maxp, ' points'
cgct         kip=pinl+maxp
cgct      endif

c  Truncate long side of interferograms if it exceed NFFT
      if(nip.ge.nfft) then
        kip=nfft-1
        if(verbose.ge.3) then
          write(6,*) 'Truncating interferogram from',nip,' to',kip
        endif
      endif
 
c  Zero-fill the uninitialized section of ac_igram up to size nfft.
c  This is necessary because the phase correction convolution is
c  performed via an FFT of size nfft.  This would not be needed
c  if we used the (slower) dot product convolution.
      call vmov(zero,0,ac_igram(kip+1),1,nfft-kip)

c  Check that the ZPD point is more than NBURST+LPCO/2 from the
c  ends of the interferogram (to avoid an array-bound error later).
      if(pinl.lt.1+nburst+lpco/2) then
        write(6,*) 'Warning from i2s_lite_siv: PINL < 1+Nburst+LPCO/2'
        write(6,*) 'ZPD is too close to end of igram: decrease LPCO'
        pinl=1+nburst+lpco/2
      endif
      if(pinl.gt.kip-nburst-lpco/2) then
        write(6,*) 'Warning from i2s_lite_siv: PINL > KIP-Nburst-LPCO/2'
        write(6,*) 'ZPD is too close to end of igram: decrease LPCO'
        pinl=kip-nburst-lpco/2
      endif

c  Best ZPD point is not necessarliy the one having the peak AC interferogram
c  amplitude. So find interferogram point closest to the maximum in symmetry.
      call bestzpd(ac_igram(pinl-nburst-lpco/2),lpco,nburst,best)
      zpdl=dble(pinl)+dble(best)
      izpd=nint(zpdl)
c      izpd=pinl   !  Old way of defining ZPD
      if(verbose.ge.3) then
        write(*,*)'Best ZPD at point: pinl,izpd= ',pinl,izpd
      endif

c      write(*,*)'izpd, lpco, nfft=',izpd, lpco, nfft
c      write(*,'(9i9)')(izpd+j,j=-4,4)
c      write(*,'(9f9.5)')(ac_igram(izpd+j),j=-4,4)

      if(ichan.eq.2) then
      call compute_odd_even_igram_shift(lgev,ac_igram(izpd-lgev/2),
     & shbar,sherr)
        izpd_si=izpd
      endif

c      if(lsemode.eq.1.or.lsemode.eq.2) then
c         if(lsemode.eq.1)stop 'lsemode=1 is not yet supported.'
c        call compute_odd_even_igram_shift(lop,
c     &  ac_igram(izpd(lsemode)-lop/2,lsemode), ! ac_igram would also need to be a 2-d array for this to work
c     &  shbar(ichan),sherr(ichan))
c      elseif(lsemode.eq.0) then
cc No LSE applied.
c        shbar(ichan)=0.0
c        sherr(ichan)=0.0
c      else
cc Eventually, we'll do different things for lsemode=3,4, but for now,
cc do nothing, and set shbar/sherr to 0.0.
c        shbar(ichan)=0.0
c        sherr(ichan)=0.0
c      endif
c      koe=(3-(-1)**izpd(lsemode))/2                        ! =1 if IZPD even; =2 if IZPD odd
c      call resample_ifg(nip,ac_igram(1:,ichan),koe-1,
c     &   shbar(ichan)*2.d0,60) ! shbar needs to be multiplied by two and converted to double for resample_ifg

      if(lsemode.eq.2) then
c Use the shbar/sherr calculated from Si.
      elseif(lsemode.eq.0) then
c No LSE applied.
       shbar=0.0
       sherr=0.0
      elseif(lsemode.eq.1) then
       stop 'lsemode=1 is not yet supported.'
      else
c Eventually, we'll do different things for lsemode=1,3,4, but for now,
c do nothing, and set shbar/sherr to 0.0.
       shbar=0.0
       sherr=0.0
      endif

      shbar_out=shbar
      sherr_out=sherr
c      write(*,*) 'ichan,shbar=',ichan,shbar,sherr,idir,izpd
c      write(59,'(7i8,2f10.5)')icall,2*(icall/4)-idir,mod((icall)/2,2),
c     & idir,izpd,(-1)**izpd,nip,shbar,sherr
c      write(59,'(8i8,2f10.6)')icall,((icall+1)/2),ichan,
c     & idir,izpd,-(-1)**(izpd+idir),(-1)**(izpd_si+idir),
c     & nip,shbar/2,sherr/2
                              
c Shift and Resample the non-ZPD-parity interferogram points by shbar/2
c The points with the same parity (odd/evenness) as IZPD are unchanged.
c SHBAR is divided by 2 because the points in BUFN have a spacing of 2
c and the SHARS subroutine treats SHBAR as the fraction of the BUFN spacing.
      nzero=nfft/2-nip/2
      nfft2=nfft/2
c      write(*,*)'nip/2,nfft/2,nzero=',nip/2,nfft/2,nzero
c     koe=(3-(-1)**izpd_si)/2                        ! =1 if IZPD even; =2 if IZPD odd
c      koe_si=(3-(-1)**izpd)/2                        ! =1 if IZPD even; =2 if IZPD odd
c      call vmov(ac_igram(koe_si),2,bufn,1,nip/2)     ! Copy non-ZPD-parity igram points to BUFN
c      call vmov(zero,0,bufn(nip/2+1),1,nzero)        ! Zero-fill remainder of BUFN
c      call shars(nfft2,bufn,-(-1)**(izpd_si+idir)*shbar/2,kflip)    ! Shift contents of BUFN by SHBAR
c      call shars(nfft2,bufn,-shbar/2,kflip)          ! Shift contents of BUFN by SHBAR
c      call vmov(bufn,1,ac_igram(koe_si),2,nip/2)     ! Merge shifted igram back into ac_igram
c      call vmov(ac_igram(koe),1,bufn,1,nip)        ! Copy igram points to BUFN
      call vmov(ac_igram,1,bufn,1,nip)        ! Copy igram points to BUFN
      if(ichan.eq.2) then  ! Si
        call resample_ifg(nip, bufn, ac_igram,
     &  1, -1.d0*dble(shbar)*(-1)**(izpd_si), 60)
      else                ! InGaAs
         call shars(1,nip,60,-(-1)**(izpd_si)*shbar,bufn,
     &   ac_igram)
      endif

c  Recompute shifts to check that they have gone to zero.
      if(ichan.eq.2) then
        call compute_odd_even_igram_shift(lgev,ac_igram(izpd-lgev/2),
     &  shbar2,sherr2)
c        write(*,*) 'ichan,shbar2=',ichan,shbar2/2,sherr2/2,idir,izpd
      endif

c      write(*,*)'izpd, lpco, nfft=',izpd, lpco, nfft

c      write(*,'(9i9)')(izpd+j,j=-4,4)
c      write(*,'(9f9.5)')(ac_igram(izpd+j),j=-4,4)

c  Copy interferogram center burst to CR in FFT-packed wrapped-around format
      do i=1,lpco
c        write(47,*)i,ac_igram(i+izpd-lpco/2)
      end do
      call vmov(ac_igram(izpd),1,cr,1,lpco/2)
      call vmov(ac_igram(izpd-lpco/2),1,cr(lpco/2+1),1,lpco/2)

c  Calculate phase correction operator
c      thresh=0.05
c      thresh=0.02   !  GCT 25-Feb-2000
      call phase_corr_oper(cr,lpco,pco_thresh,phasepath)

c  Perform in-place convolution of raw interferogram in RBUF with CR
c    Before convolution, igram extends from 1 to kip
c    After, usable igram extends from 1+lpco/2 to kip-lpco/2+1 (KIP-LPCO+1 pts)
      call convol(cr,lpco,ac_igram,nfft)
ct      call symmetry(ac_igram(izpd),zsym,symmp,lpco)  ! corrected igram symm
ct      write(*,*) ' Corrected ZPD: Index, Value, Symmetry  =',
ct     $ izpd,ac_igram(izpd),zsym
ct      zpdv=ac_igram(izpd)
ct      totp=kip-izpd-lpco/2+1

c  Move corrected interferogram so that ZPD is at 1 instead of PINL
c  overwriting short side of interferogram which is lost.
      call vmov(ac_igram(izpd),1,ac_igram(1),1,kip-lpco/2+1-izpd) ! dg
c  Interferogram now extends from 1 to KIP-LPCO/2+1-IZPD.

cgctc  Multiply phase-corrected igram by a windowing function.
cgct      me=1.0
cgct      do i=1,kip-lpco/2+1-izpd
cgct       ac_igram(i)=ac_igram(i)*(1+(me-1)*float(i-1)/(kip-lpco/2+1-izpd))
cgct      end do

      zpa=dble(ac_igram(1))  !  ZPD interferogram amplitude (phase-corrected)

c  Zero-fill unused portions of RBUF beyond long side of interferogram,
c  including the NFFT+1'st point which otherwise would not subsequently
c  get initialized when the interferogram is unfolded.
      call vmov(zero,0,ac_igram(kip-lpco/2+2-izpd),1,
     $  nfft-kip+lpco/2+izpd)  ! dg

c  Unfold long side of interferogram, together with appended DC-filled points
c      call vmov(ac_igram(2),1,ac_igram(2*nfft),-1,nfft-1)
      j=2*nfft
      do i=2,nfft
         ac_igram(j)=ac_igram(i)
         j=j-1
      end do

c  Real-to-complex FFT of double-sided, perfectly symmetrical interferogram
      call ffak(ac_igram,2*nfft)

c  Perform compaction (i.e. remove imaginary values which are all zero anyway)
c  and normalize spectral values by dividing by 2*NFFT (multiply by 0.5/NFFT)
c     call vmul(ac_igram,2,0.5/nfft,0,ac_igram,1,nfft)
      call vmul(ac_igram,2,0.125,0,ac_igram,1,nfft)
c  Above line changed to match OPUS values (we were NFFT/4 too small)

      nip=nfft   ! Return spectrum size

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
c     RBUF      the raw interferogram
c     NP        number of points over which to compute symmetry
c     NBURST    number of interferogram points to be tested
c
c  OUTPUTS:
c     BEST      fractional interferogram point having highest symmetry
c
c  Note that the symmetry is evaluated every half point, and then the 3
c  half points bracketing the maximum symmetry are fitted by a parabola
c  in order to determine the fractional location of peak symmetry.
c  The best point to use would usually be NINT(BEST).

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
c==================================================================
      subroutine symmetry(ac_igram,symmi,symmp,np)
c
c  Computes the symmetry of an interferogram about the point RBUF(0)
c  and also about the mid-point between RBUF(0) and RBUF(1).
c
c  INPUTS:
c     NP                 number of points over which symmetry is to be evaluated
c     RBUF(-np/2:np/2)   section of raw interferogram
c
c  OUTPUTS:
c     SYMMI              the interferogram symmetry about RBUF(0.0)
c     SYMMP              the interferogram symmetry about RBUF(0.5)
c
c  Symmetry is evaluated by comparing the sums and differences of points
c  which are symmetrical about RBUF(0). If the interferogram is perfectly
c  symmetrical, then  SADELI = SUM{ABS[RBUF(-i)-RBUF(+i)]} = 0, whereas
c  if igram is antisymmetrical then  SASUMI = SUM{ABS[RBUF(-i)+RBUF(+i)]} = 0.
c  Therefore, symmetry can be defined by (SASUMI-SADELI)/(SASUMI-SADELI), which
c  will vary from -1 for perfect anti-symmetry, to +1 for perfect symmetry.
c
c  Similary, the symmetry can also be defined about the half point RBUF(0.5)
c  such that SADELP = SUM{ABS[RBUF(-i+1)-RBUF(+i)]} and
c            SASUMP = SUM{ABS[RBUF(-i+1)+RBUF(+i)]}

      implicit none

      integer*4 np,j
      real*4 ac_igram(-(np/2):+np/2),
     & sasumi,sadeli,sasump,sadelp,symmi,symmp,q,x,ww,pi

      pi = real(4.d0*datan(1.d0))
      sasumi=0.0
      sadeli=0.0
      sasump=0.0
      sadelp=0.0
      q=pi/float(np)
      do j=1,np/2
        x=q*float(j)
        ww=(5.0*cos(x)+cos(3.0*x))/6.0
        sasumi=sasumi+ww*abs(ac_igram(-j)+ac_igram(j))
        sadeli=sadeli+ww*abs(ac_igram(-j)-ac_igram(j))
        sasump=sasump+ww*abs(ac_igram(-j+1)+ac_igram(j))
        sadelp=sadelp+ww*abs(ac_igram(-j+1)-ac_igram(j))
      enddo
      symmi=(sasumi-sadeli)/(sasumi+sadeli)
      symmp=(sasump-sadelp)/(sasump+sadelp)
      return
      end
c==================================================================
      subroutine phase_corr_oper(cr,lpco,pco_thresh,phasepath)
c  Takes the interferogram center burst (CR) and calculates an operator which
c  when convolved with the raw interferogram, should make it symmetrical.
c 
c  First apodizes the raw interferogram center burst in array CR(LPCO).
c  Then computes the complex low resolution spectrum (a+i.b) where i=sqrt(-1).
c
c  The complex low-resolution spectrum contains LPCO/2+1 real values (a) and
c  LPCO/2-1 imaginary values (b). The first and last b values, representing
c  the DC and Nyquist frequencies, are zero by definition because the antisymetric
c  component of the interferogram must have DC and Nyquist values of zero.
c  These values are omitted from the returned spectrum and the first imaginary
c  element of theretunred spectrum actually contains the real Nyquist amplitude,
c  not the DC imaginary.
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
c             LPCO           Length of Phase Correction Operator
c             CR(LPCO)       Array containing interferogram center burst
c             pco_thresh     The fractional intensiy threshold for use of phase information
c             phasepath      Where to write the phase curve
c
c OUTPUTS:
c             CR(LPCO)       Array containing phase correction operator

      implicit none

      integer*4 mpco,lpco,i,lastoki,lnbc
      parameter (mpco=64*1024)
      real*4 cr(lpco),amplit(1+mpco/2),phase(1+mpco/2),area,rarea,
     & ampmax,pythag
c      real*4 crodd(mpco/2),creven(mpco/2)
      real*8 pco_thresh,flag,freq
      character phasepath*(99)
c
      if(lpco.gt.mpco) stop 'Error in phase_corr_oper: LPCO > MPCO'

      lastoki=0
      flag=1.0d0
      if(cr(1).lt.0.0) flag=-1.d0
      call apodize(cr,lpco)

c      call vmov(cr(1),2,crodd,1,lpco/2)
c      call vmov(cr(2),2,creven,1,lpco/2)
c
c  Perform small real-to-complex FFT of apodized interferogram center burst
      call ffak ( cr, lpco )

c      call ffak ( crodd, lpco/2 )
c      call ffak ( creven, lpco/2 )
cc  Write out phase and amplitude curves
c      lp=lnbc(phasepath)
c      if(lp.gt.0) then
c         open(19,file=phasepath(:lp)//'.odd')
c         write(19,*)2,3
c         write(19,*)' Freq  Amplit  Phase '
c         do i=1,lpco/4
c            freq=15798*float(i-1)/(lpco/2)
c            amplit(i)=sqrt(crodd(2*i-1)**2+crodd(2*i)**2)
c            ww=flag/amplit(i)
c            phase(i)=datan2 ((-ww)*crodd(2*i), ww*crodd(2*i-1) )
c            write(19,'(f9.1,3f9.4)') freq,amplit(i),phase(i)
c         end do
c         close(19)
c         open(19,file=phasepath(:lp)//'.even')
c         write(19,*)2,3
c         write(19,*)' Freq  Amplit  Phase  '
c         do i=1,lpco/4
c            freq=15798*float(i-1)/(lpco/2)
c            amplit(i)=sqrt(creven(2*i-1)**2+creven(2*i)**2)
c            ww=flag/amplit(i)
c            phase(i)=datan2 ((-ww)*creven(2*i), ww*creven(2*i-1) )
c            write(19,'(f9.1,3f9.4)') freq,amplit(i),phase(i)
c         end do
c         close(19)
c      endif
c
c  Compute low resolution power and phase spectra,
      ampmax=0.0d0
c      tot=0.0d0
c      toty=0.0d0
      cr(1)=abs(cr(1))     ! Set DC magnitude 
      cr(2)=0.0            ! Set DC phase 
      amplit(1)=cr(1)      ! DC
      phase(1)=0.0         ! DC
      do 40 i=2,lpco/2
         freq=15798*float(i-1)/(lpco/2)
c        mag=sqrt(cr(2*i-1)**2+cr(2*i)**2)
        amplit(i)=pythag(cr(2*i-1),cr(2*i))
c       write(*,*)'pythag,orig=',amplit(i),
c    & sqrt(cr(2*i-1)**2+cr(2*i)**2)
c        amplit(i)=sqrt(cr(2*i-1)**2+cr(2*i)**2)
        if(amplit(i).gt.ampmax) ampmax=amplit(i)
        if(amplit(i).gt.pco_thresh*ampmax) lastoki=i
c        ww=flag/amplit(i)
c        theta=datan2 ((-ww)*cr(2*i), ww*cr(2*i-1) )
c        phase(i)=datan2 ((-ww)*cr(2*i), ww*cr(2*i-1) )
        phase(i)=datan2 ((-flag)*cr(2*i), flag*cr(2*i-1) )
c        cr(2*i-1)=sngl(mag)
c        cr(2*i)=sngl(theta)
        cr(2*i-1)=amplit(i)
        cr(2*i)=phase(i)
c        toty=toty+theta*mag
c        tot=tot+mag
40    continue
c
c  Interpolate phase across regions will little power (mag < pco_thresh*ampmax)
c      call vintrp(cr,2,lpco/2,lastoki,1)
      if(lastoki.eq.0) then
        write(*,*)'Warning: lastoki being set to lpco/2. This is '//
     & 'probably a bad spectrum.'
        lastoki=lpco/2
      endif
      call vintrp(phase,1,lpco/2,lastoki,1)
      lastoki=1
      do 41 i=lpco/32,lpco/2
      if(amplit(i).gt.pco_thresh*ampmax) then
         if(i.ne.lastoki+1) call vintrp(phase,1,lpco/2,lastoki,i)
         lastoki=i
      endif
41    continue
c
c  Write out phase and amplitude curves
      if(lnbc(phasepath).gt.0) then
         open(19,file=phasepath)
         write(19,*)2,4
         write(19,*)' Freq  Amplit  Phase1  Phase2 '
         do i=1,lpco/2
            freq=15798*float(i-1)/(lpco/2)
c           write(19,'(f9.1,f9.4,f11.2)')2*15798*float(i-1)/lpco,
           write(19,'(f9.1,3f9.4)') freq,amplit(i),cr(2*i),phase(i)
         end do
         close(19)
      endif
c
c  Compute complex exponential of the adjusted phase
      do i=1,lpco/2
         cr(2*i-1)=cos(phase(i))
         cr(2*i)=sin(phase(i))
      end do
      cr(2)=1.0   ! insert the Nyquist frequency at cr(2)
c
c  take the inverse transform (complex-to-real)
      call ffsk ( cr, lpco )
c
c  Apodize phase correction operator
      call apodize( cr, lpco )
c
c  normalize the phase-correcting function to unit area
      call vdot(cr,1,1.0,0,area,lpco)
      rarea=sngl(flag)/abs(area)
      call vmul(cr,1,rarea,0,cr,1,lpco)
c      sum=0.
c      do i=1,lpco
c        sum=sum+cr(i)
c      enddo
c      t=flag/abs(sum)
c      do i=1,lpco
c        cr(i)=t*cr(i)
c      enddo
      return
      end
c
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

      subroutine apodize(y,lpco)
c   Performs in-place apodization of array Y according to the function 
c   ( 5*cos(i*pi/LPCO) + cos(3*i*pi/lpco) ) / 6
c   Assumes that Y is packed ready for FFT, i.e. with +ve part of operator
c   in Y(1:LPCO/2) and -ve half of operator in Y(1+LPCO/2:LPCO)

      implicit none

      integer*4 lpco,i
      real*4 y(lpco),x,pi,q,ww
c
      pi = real(4.d0*datan(1.d0))
      q=pi/float(lpco)
      do i=1,lpco/2-1
        x=q*float(i)
        ww=(5.0*cos(x)+cos(3.0*x))/6.0
        y(1+i)=y(1+i)*ww
        y(lpco+1-i)=y(lpco+1-i)*ww
      enddo
      y(1+lpco/2)=0.0
c  Note that Y(1) is unchanged
      return
      end
c==================================================================
      subroutine convol(cr,lpco,ac_igram,nfft)
c  Performs a fast convolution of the raw interferogram in RBUF with the phase
c  correction operator CR.
c  The resulting phase corrected interferogram is written in-place in RBUF,
c  with its ZPD in the same location as that of the raw interferogram.
c
c  INPUTS:
c     LPCO           Length of Phase Correction Operator
c     NFFT           Size of FFT
c     CR(LPCO)       Phase correction operator
c     RBUF(NFFT)     Raw interferogram
c
c  OUTPUTS:
c     CR(LPCO)                    Unchanged
C     RBUF(1:LPCO/2-1)            Garbage (spectral contamination)
C     RBUF(LPCO/2:NFFT-LPCO/2)    Corrected interferogram (NFFT-LPCO+1 points)
C     RBUF(NFFT-LPCO/2+1,2*NFFT)  Garbage (spectral contamination)
c
c  Note that the declared dimension of RBUF is 2*NFFT because the upper half
c  of this array is used as workspace (to FFT the phase correction operator).
c  Also note that on exit, even the lower half of RBUF will contain garbage
c  outside the range LPCO/2 to NFFT-LPCO/2 due to wrap-around effects.
c
c  Note that there are two possible ways of doing this convolution:
c    1) In the time domain, which is simple but slow
c    2) In the spectral domain, which is messy (because it requires forward &
c       reverse FFT's), but fast (it requires only a complex multiplication).
c  This subroutine uses the latter method, but also contains the equivalent
c  code (commented) for the time domain, which has been verified to give
c  identical results.

      implicit none

      integer*4 lpco,nfft,i,j
      real*4 ac_igram(2*nfft),cr(lpco),xr
c
c================================================================
cc  Perform time-domain (slow) convolution.
cc  Note that reversal of CR places center point at LPCO/2 rather than LPCO/2+1
c      do i=1,nfft-lpco+1
c         call vdot(ac_igram(i),1,cr(lpco/2),-1,x1,lpco/2)
c         call vdot(ac_igram(i+lpco/2),1,cr(lpco),-1,x2,lpco/2)
c         ac_igram(i)=x1+x2
c      enddo
cc  Move corrected igram to line up with raw igram
c      call vmov(ac_igram(nfft-lpco+1),-1,ac_igram(nfft-lpco/2),-1,nfft-lpco+1)
cc  Be aware that this leaves garbage in RBUF locations 1 to LPCO/2
cc================================================================
c  Move CR to the ends of upper half of RBUF and set middle to zero
c      call vmov(cr,1,ac_igram(nfft+1),1,lpco/2)
c      call vmov(cr(lpco),-1,ac_igram(2*nfft),-1,lpco/2)
c      call vmov(0.0,0,ac_igram(nfft+1+lpco/2),1,nfft-lpco)
      do j=1,lpco/2
         ac_igram(nfft+j)=cr(j)
         ac_igram(2*nfft+1-j)=cr(lpco+1-j)
      end do
      do j=nfft+lpco/2+1,2*nfft-lpco/2
         ac_igram(j)=0.0
      end do

c
c  Perform a full size FFT in order to interpolate complex phase spectrum.
      call ffak(ac_igram(nfft+1),nfft)
c
c  Perform real-to-complex transform on raw interferogram in lower half of RBUF.
      call ffak(ac_igram,nfft)
c
c  Multiply resulting complex spectrum by interpolated complex phase spectrum
c  in upper half of RBUF
      ac_igram(1)=ac_igram(1)*ac_igram(nfft+1)   !  No imaginary DC term
      ac_igram(2)=ac_igram(2)*ac_igram(nfft+2)   !  No imaginary Nyquist term
      do i=4,nfft,2
        j=nfft+i
        xr=      ac_igram(i-1)*ac_igram(j-1)-ac_igram(i)*ac_igram(j)
        ac_igram(i)= ac_igram(i-1)*ac_igram(j)+ac_igram(j-1)*ac_igram(i)
        ac_igram(i-1)=xr
c        cbuf(i/2)=cbuf(i/2)*cbuf(nfft/2+i)  ! equivalent to 4 previous lines
      enddo
c
c  Inverse FFT complex spectrum to yield corrected igram in lower half of RBUF.
      call ffsk(ac_igram,nfft)
      return
      end
