      subroutine abscoj(runlab,nlev_mav,t,p,d,nspeci,targmol,
     & wmf,vpf,
     & tll_file,parfile,ifcsp,grid,ncp,vac,vvoigt,calc_lm,calc_nv)
c
c  Performs line-by-line calculations of the volume absorption coefficients,
c  which are stored in a 3-D array (frequency, altitude, target gas).
c  First it reads "molparam.dat" and pre-computes the partition functions.
c  Then it reads the linelists and does the line-by-line calculations.
c
c  INPUT Parameters:
c    nlev_mav             I*4   Number of atmospheric levels.
c    T(nlev_mav)          R*4   Temperature Profile (Kelvin)
c    P(nlev_mav)          R*4   Pressure Profile (Atmospheres)
c    D(nlev_mav)          R*4   Number Density Profile (molec.cm-3)
c    nspeci               I*4   Number of different species/isotopomers in PARFILE
c    targmol(nspeci)      I*4   Group to which absorption of each specie belongs.
c    wmf(nspeci,nlev_mav) R*4   Array of Wet Mole Fraction (WMF) profiles for each gas
c    vpf(nspeci,nlev_mav) R*4   Workspace Array (to store partition functions)
c    tll_file             C**   Spectral linelists to be read.
c    parfile              C**   File of molecular parameters (e.g. molparam.dat).
c    fzero                R*8   Frequency immediately prior to the first grid point.
c    grid                 R*8   Primitive grid point spacing (cm-1).
c    ncp                  I*4   Number of primitive spectral grid points.
c    vvoigt(ncp)          R*4   Work Array used to store Voigt lineshape
c    calc_lm              L*1   Logical saying whether to do LM or not
C    calc_nv              L*1   Logical saying whether use Non Voigt  
c
c  OUTPUT Parameter:
c    vac(ncp,nlev,*)      R*4   3-D matrix of absorption coefficients (cm-1)
c
c
c  VAC contains the Absorption Coefficient spectrum for each level for
c  each target gas in units of cm-1. It is the sum of the line absorptions
c  in units of cm-1/(molec.cm-1), multiplied by the gas Wet Mole Fraction
c  (wmf) and total number density.  Its dot product with the vector of slant
c  paths in cm yields the optical thickness (dimensionless).
c  Absorption coefficients are stored at frequencies vj=fzero+j.grid  j=1,ncp
c  at nlev different levels, and for different groups of species.
c  Note that in the main program VAC is stored as a 1-D array and so we
c  don't have to worry here about its logical versus physical dimensions.
c
c  Line shape is truncated when the monochromatic optical depth falls below 
c  To speed up calculation of weak or narrow lines, the line-shape is
c  truncated. To prevent this causing discernable discontinuities in the
c  calculated transmittance spectrum, the truncation is delayed until
c  the monochromatic Optical Depth fall below OD' (= 1E-04).
c  In the far line wings, the lineshape is Lorentzian, and so the
c  monochromatic OD can be approximated by
c      OD(x) = EW. w/Pi/[w^2+x^2]
c  where x is the distance (cm-1) from line center, w is the Lorentz width,
c  and EW is the equivalent width of the line.
c      x = w.Sqrt[EW/Pi/OD'/w-1]
c  So OD(x')  = OD' when  x'= w.Sqrt[EW/Pi/OD'/w - 1]
c  To prevent Sqrt[-ve] when EW small...........
c
c  Another consideration is not to truncate lines that still have
c  sinfificant integrated OD across the window, even though the
c  monochromatic OD may have already fallen below OD'. For example,
c  there may be dozens of lines, some located outside the window,
c  whose individual absorptions may be less than OD', but whose
c  collective contribution to the transmittance is large. We want
c  to retain these lines until their remaining equivalent width
c  (integrated OD) falls below EW".  The integral of the pure Lorentz
c  lineshape is EW/Pi.tan-1(x/w).  In the limit that x is close to
c  infinity, this can be approximated by  EW.(0.5-w/Pi/x), so the
c  line area between x=x" and x=Inf is EW.w/Pi/x".  So we keep on
c  using a line until the remaining area is less than EW", which
c  is until x" = EW.w/Pi/EW"
c
c  So the line-shape integration continues until both of the
c  following criteria are fulfilled:
c      x' > w.Sqrt[EW/Pi/OD'/w - 1]
c  or  x" > EW.w/Pi/EW"
c  OD' is the monochromatic OD threshold (OD'=1E-4)
c  EW" is the integrated OD threshold (EW"=5E-5)
c
c  Note that the former criterion varies as Sqrt(EW.w) whereas the
c  latter varies as EW.w.  So for strong, broad lines the latter
c  criterion will tend to dictate the limits of the line-shape
c  integration, whereas for weaker narrower lines the former
c  criterion will delay truncation
c
c  For a strong broad H2O line with w=0.1 cm-1 and EW=0.1 cm-1
c  then the two criteria yield
c      x' > 0.1.Sqrt[1000/Pi - 1]   = 1.8 cm-1
c      x" >  200/Pi                 = 64 cm-1
c
c  For a strong O3 line with w=0.005 and EW=0.1 the two criteria yield
c      x' > 0.005.Sqrt[1000.Pi - 1] = 0.9 cm-1
c  or  x" > 10/Pi =                 = 3.2 cm-1
c    
c  For a weak tropospheric H2O line with w=0.1 cm-1 and EW=0.01 cm-1
c  then the two criteria yield
c      x' > 0.1.Sqrt[1000/Pi - 1]   = 1.8 cm-1
c  or  x" > 20/Pi                   = 6.4 cm-1

c  For a weak stratospheric O3 line with w=0.005 cm-1 and EW=0.001
c      x' > 0.005.Sqrt[2000/Pi - 1] = 0.125 cm-1
c  or  x" > 0.1/Pi                  = 0.032 cm-1

      implicit none

      real*4  spi         !PI in single precision
      parameter(spi=3.14159265)

      include "../gfit/ggg_int_params.f"

      integer*4
     &   nlhead,ncol_iso,istat,ciso2kiso,
     &   kline1,kline2,jline,nlev_mav,
     &   ncp,lnbc,lloc,llf,lr,j,
     &   atmnv_flag,
     &   irec,lq,
     &   jlo,jhi,
     &   klo,khi,
     &   ifcsp,
c     &   jtg, jcp,
     &   nltnu,       ! total number of lines exceeding TNULST
     &   lunr_ll,      ! parfile & linelists
     &   lunr_iso,    ! isotopologs
     &   lunw_rpt,    ! 'abscof.rpt'
c     &   lunw_vac,     ! 'vac_ilev01.out'
     &   lunr_lm_bin, ! reading binary LM files
     &   kv1, kv2,    ! limits of line profile calculation
     &   nv,jv,       ! number of Voigt point = kv2 - kv1 +1
     &   jspeci,      ! Indexes the different isotopomers
     &   ispeci,      ! Indexes the different isotopomers
     &   nspeci,      ! Number of different species/isotopomers in PARFILE
     &   nvmode,      ! Number of different vibrational modes (ie 3N-6) of each gas
     &   ngas,kgas,jgas,kwas,  ! Number of different gases found in MOLPAR
     &   idot,
     &   mlf,nlf,     ! maximum & actual number of linefiles
     &   posnall,jtarg,mm,mmax,reclen,nlines,kiso,i,k,
     &   lmax,lfile,ilev,idum,lrt,lp
c     &   np           ! # of lines in relaxation matrix needed by full LM

      integer*8 fsib,file_size_in_bytes
c
      parameter (mlf=99,lunr_ll=21,lunr_iso=22,lunw_rpt=23,
c     & lunw_vac=24,
     & lunr_lm_bin=46)
      integer*4
     &   targmol(nspeci),
     &   specindex(mgas+1), ! specindex(kgas)+kiso  is the specie number/identifier
     &   gasindex(nspeci),  ! gasindex(ispeci) is the gas number
     &   maxlev(nspeci),    ! level at which each species has its max density
     &   dgen(mvmode),      ! degeneracy of vibrational modes
     &   icode,             ! Isotopolog code
     &   molewt(nspeci)     ! Molecular Weight
c
      real*8
     &   fzero,         ! frequency immediately prior to first grid point
     &   fmax,          ! fzero+ncp*grid
     &   fcen,          ! fzero+ncp*grid/2
     &   stren,         ! line strength
     &   hw,            ! half-width (cm-1) of precomputed interval
     &   nuoff,         ! search linelist NUOFF cm-1 beyond ends of VAC array 
     &   nu1,nu2,       ! start & end frequencies (cm-1) for linelist search
     &   grid,          ! primitive point spacing (cm-1)
     &   freq,          ! line center frequency (cm-1)
     &   delpl(nspeci), ! frequency spacing of consecutive pseudo-lines
     &   fprev(nspeci), ! line center frequency of previous line of JSPECI
     &   vcent,         ! line center position in primitive grid points
     &   cen,           ! central frequency of window (cm-1)
     &   x1,            ! start frequency  (for Voigt evaluation)
     &   godw,          ! GRID / Doppler-width
     &   flinwid,       ! half-width of computed lineshape
     &   Yair1,         ! Line-mixing air parameter 1
     &   Yair2,         ! Line-mixing air parameter 2
     &   Yair3,         ! Line-mixing air parameter 3
     &   Yair4,         ! Line-mixing air parameter 4
     &   Yself1,        ! Line-mixing self parameter 1
     &   Yself2,        ! Line-mixing self parameter 2
     &   Yself3,        ! Line-mixing self parameter 3
     &   Yself4,        ! Line-mixing self parameter 4
     &   Yh2o1,         ! Line-mixing H2O parameter 1
     &   Yh2o2,         ! Line-mixing H2O parameter 2
     &   Yh2o3,         ! Line-mixing H2O parameter 3
     &   Yh2o4,         ! Line-mixing H2O parameter 4
     &   YT,            ! First order line-mixing parameter
     &   sg,            ! First grid point in cm-1
     &   HWT,           ! width to pass to qSDV
     &   SDHWT,         ! SDwidth to pass to qSDV
     &   SHIFT,         ! p-shift to pass to qSDV
     &   SDSHIFT        ! SD p-shift to pass to qSDV 
c
      real*4
     &   AbsY,
     &   AbsV,AbsW,
     &   ewvb(nspeci),    ! Extra water vapor broadening
     &   vvoigt(ncp),   ! Array of Voigt profile
     &   y,             ! Lorentz-width/Doppler-width (used by VOIGT(x,y))
     &   tnu,tnulst,    ! line-center absorbtance & threshold
     &   linwid,        ! approximate line width (cm-1)
     &   del,           ! line center distance beyond edge of VAC interval
     &   srpi,          ! SQRT(3.14159265...)
     &   fia(nspeci),   ! Fractional Isotopolog Abundance
     &   delta,lnfrd, ! Fractionation at tropopause and its gradient
     &   eprime,        ! ground state energy (cm-1)
     &   eptf,          ! eprime*tfac
     &   abhw, abhw_h2o,! air-broadened half-width
     &   sbhw,          ! self-broadened half width
     &   pbhw,          ! pressure-broadened half width
     &   sdwp,          ! speed-dependent width parameter
     &   sdshiftp,      ! speed-dependent shift parameter
     &   dn,            ! dicke narrowing parameter
     &   corr,          ! correlation for collisional narrowing
     &   cgs,cgsmax,    ! concentration of ground-state molecules
     &   pshift,        ! line-center frequency pressure shift (cm-1/atm)
     &   spshift,       ! line-center frequency self press shift (cm-1/atm)
     &   trat(mlev),    ! trat(k)=296./t(k)
     &   tcurv(mlev),   ! tcurv(k)=(1.-t(k)/280.)*(t(k)/220.-1.)/0.01461)
     &   tfac(mlev),    ! tfac(k)=conx*(trat(k)-1.0)
     &   stimem(mlev),  ! stimem(k)=(1.-exp(sf*trat(k)))/(1.-exp(sf))
     &   atc(nspeci),   ! Additional Temperature Dependence (PLL)
     &   tdrpf(nspeci), ! T-Dependence of Rotational Partition Function
     &   tdpbhw,        ! T-Dependence of Pressure-Broadened Half-Width
     &   tdpshift,      ! T-Dependence pressure shift (units?)
     &   tdsbhw,        ! T-Dependence of Self Pressure-Broadened Half-Width
     &   tdspshift,     ! T-Dependence self pressure shift (units?)
     &   vibfrq(mvmode),! Array of vibrational frequencies
     &   conx,          ! conx=-1.43881d0/296  ( = hc/kT;  T=296)
     &   dopphwem,      ! Exact HWEM Doppler width
     &   sublor,        ! line width after which sub-Lorentzian wings kicks in
     &   alor,          ! actual Lorentz width (cm-1). Pressure * PBHW
     &   frac,
     &   vibpf,         ! Vibrational Partition Function.
     &   sf,            ! scratch variable used in computation of STIMEM
     &   sxcgs,         ! scratch variable used to save STREN*CGS
     &   sxcgsorpidw,   ! scratch variable multiplying VOIGT in innermost loop
     &   t296,          ! =296.0
     &   v296,          ! 296 K vibrational partition function
     &   wmf_tot(mgas,nlev_mav), ! Total FIA-weighted WMF for each gas
     &   wmf(nspeci,nlev_mav),   ! Wet Mole Fraction
     &   vpf(nspeci,nlev_mav),vmax,
     &   vac(ncp,nlev_mav,0:*),
     &   t(nlev_mav),p(nlev_mav),d(nlev_mav)
      parameter (t296=296.,conx=-(1.43881/t296))

      character
     & iso_fmt*(80),           ! Format of isotopologs.dat
     & runlab*(*),            ! spectrum name
     & lmfile*120,            ! file containing LM absorption coefficients
c     &   vacfile*48,
c     &   rcsmax*40,
     &   fmt_rpt*80,
     &   llformat*49,         ! FORMAT statement for linelist
     &   llformat2*155,       ! FORMAT statement for the atmnv linelist
     &   ciso*1,
     &   quantum*94,          ! transition quantum numbers
     &   gasname*8,           ! Gas name
     &   linfil(mlf)*(mfilepath),  ! Name of currently open linelist
     &   parfile*(*),         ! Name of MOLPAR file (usually 'isotopologs.dat")
     &   speci_id(nspeci)*20, ! chemical formulae of the species/isotopomers
     &   tll_file*(*),        ! name of file containing names of linelists
     &   sfile_path*(mpath+18),! Path to file needed for Full LM
     &   sfile*14,            !Full LM file name
     &   gggdir*(mpath),dl*1  !GGG directory path
c
      logical scia,fcia,hitflag,linfil_exist,lmfile_exist,calc_lm,
     &        calc_nv
c====================================================================
      idum=mauxcol  ! Avoid compiler warning (unused parameter)
      idum=mcolvav  ! Avoid compiler warning (unused parameter)
      idum=mrow_qc  ! Avoid compiler warning (unused parameter)
      idum=mspeci   ! Avoid compiler warning (unused parameter)
      idum=ncell    ! Avoid compiler warning (unused parameter)
      idum=nchar    ! Avoid compiler warning (unused parameter)

      fmt_rpt='(f10.3,2i3,i4,1pe10.2,0pf8.1,f5.3,2x,a18,f7.1,i4,i8,i4)'

c Extra water vapor broadening (compared with air broadeniong).
c      do i=1,mgas
c         ewvb(i)=0.4 ! default value
c      end do
cc For H2O HDO, and D2O this is already included in the SBHW
c      ewvb(1) = 0.0
c      ewvb(49) = 0.0
c      ewvb(71) = 0.0
c      ewvb(2) = 0.7
c      ewvb(4) = 0.6
c      ewvb(5) = 0.7
c      ewvb(6) = 0.36
c      ewvb(7) = 0.40
c      ewvb(28) = 2.6

      mmax=0
      lq=0   ! Avoid compiler warning
      srpi=sqrt(spi)
c      write(*,*)'abscoh: t/p/d = ',t(1),p(1),d(1)
      if(nlev_mav.gt.mlev) then
         write(*,*)'nlev,mlev=',nlev_mav,mlev
         stop 'ABSCO: Increase parameter MLEV'
      endif
c      if(nspeci.gt.mspeci) stop 'ABOCOF: Increase parameter MSPECI'
      fzero=grid*(ifcsp-1)
      fmax=fzero+ncp*grid
      fcen=fzero+ncp*grid/2 
      open(lunw_rpt,file='abscof.rpt',status='unknown')
      write(lunw_rpt,43)' Pre-computing VACs from ',
     & fzero+grid,'  to ',fmax,' cm-1 :',ncp,
     & ' grid points at a primitive point spacing of ',grid, ' cm-1'
43    format(a,f12.6,a,f12.6,a,/,i6,a,f13.10,a)
c
c   SET PARAMETERS FOR TRANSMISSION CALCULATION
      hw=grid*(ncp-1)/2
      cen=fzero+grid+hw
      frac=5.0E-06   !  capture 99.9995% of line absorption
      nuoff=10.D0*dsqrt(2+dble(grid)*ncp*dble(vmax(p,1,nlev_mav)))
      if(nuoff.gt.700) nuoff=700.0  ! Venus
c      write(*,*)'nuoff=',nuoff,grid*ncp
      
      nu1=fzero+grid-nuoff
      nu2=fmax+nuoff 
c
c  PRE-CALCULATE TRAT, TFAC, & STIMULATED EMISSION TERM FOR EACH LAYER
c  Assume that the central window wavenumber CEN is close to the line.
      sf=conx*sngl(cen)
      do k=1,nlev_mav
         if(t(k).lt.50.) write(*,*)'abscoj: warning: T < 50K',t(k)
         trat(k)=296./t(k)
         tfac(k)=conx*(trat(k)-1.0)
         stimem(k)=(1.-exp(sf*trat(k)))/(1.-exp(sf))
         tcurv(k)=(1.-t(k)/280.)*(t(k)/220.-1.)/0.01461
      end do
c-------------------------------------------------------------------
c  CALCULATE THE VIBRATION PARTITION FUNCTIONS (VPF)
      open(unit=lunr_iso,file=parfile(:lnbc(parfile)),status='old')
      read(lunr_iso,*) nlhead,ncol_iso
      read(lunr_iso,'(7x,a)') iso_fmt
      do j=3,nlhead
         read(lunr_iso,*)               ! Column headers
      end do
      do jspeci=1,nspeci
         call read_isotopolog(lunr_iso,iso_fmt,kgas,kiso,gasname,
     &   speci_id(jspeci),icode,fia(jspeci),delta,lnfrd,molewt(jspeci),
     &   ewvb(jspeci),atc(jspeci),tdrpf(jspeci),vibfrq,dgen,nvmode,
     &   mvmode,istat)
c         write(*,*) jspeci,kgas,kiso,gasname,speci_id(jspeci),icode
         if(istat.ne.0) exit
         if(kgas.lt.0) stop 'KGAS<0'
         if(kiso.lt.0) stop 'KISO<0'
         gasindex(jspeci)=kgas
         specindex(kgas)=jspeci-kiso
         if( specindex(kgas).lt.0) then
            write(*,*) jspeci,kgas,kiso, specindex(kgas)
            stop 'specindex < 0'
         endif

c
c  CALCULATE VPF AT 296K
         v296=vibpf(t296,vibfrq,dgen,nvmode)  ! VPF at 296 K
c
c  CALCULATE VPFs AT nlev OTHER TEMPERATURES AND DIVIDE BY v296
c         write(*,*) ' 1: kgas,kiso,jspeci,gasindex(jspeci) = ',
c     &   kgas,kiso,jspeci,gasindex(jspeci)
         do k=1,nlev_mav
            vpf(jspeci,k)=vibpf(t(k),vibfrq,dgen,nvmode)/v296  ! VPF(T)/VPF(296)

c  Fudge T-dependence of N2O5, CF4, CFC-11, CFC-12, HCFC-22, CFH2CF3.
c  The factor 0.01461 = 1/((1-250/280)*(250/220)) built into tcurv
c  normalizes the  ATC value to be the correction at 250K.
            vpf(jspeci,k)=vpf(jspeci,k)*(1.-atc(jspeci)*tcurv(k))
         end do

      end do  ! jspeci=1,nspeci
      close(lunr_iso)
      ngas=kgas
      specindex(ngas+1)=jspeci
      if(ngas.gt.mgas) STOP ' ABSCOF: Increase parameter MGAS'

c  Compute total DMF for each gas (used for self-broadening).
c  This is done by adding the true WMFs (FIA*WMF) for each isotopolog
      do ilev=1,nlev_mav   ! Loop over Levels
         jgas=1
         jspeci=1
         wmf_tot(jgas,ilev)=fia(jspeci)*abs(wmf(jspeci,ilev))
         do jspeci=2,nspeci   ! Loop over Gases
c  Sum over isotopologs of each gas
            if(jgas.eq.gasindex(jspeci)) then
c            write(*,*)jspeci,jgas,gasindex(jspeci)
               wmf_tot(jgas,ilev)=wmf_tot(jgas,ilev)+
     &         fia(jspeci)*abs(wmf(jspeci,ilev))
            else
               jgas=jgas+1
               wmf_tot(jgas,ilev)=fia(jspeci)*abs(wmf(jspeci,ilev))
            end if
c            if(ilev.eq.3) write(*,*) jspeci,jgas,fia(jspeci),
c     &      wmf(jspeci,ilev),wmf_tot(jgas,ilev)
         end do
c  Water vapour is a special case, H2O, HDO & D2O have to be added.
         wmf_tot(1,ilev)=
     &   wmf_tot(1,ilev)+wmf_tot(49,ilev)+wmf_tot(71,ilev)     ! H2O
c     &   9.653449e-08*wmf_tot(71,ilev)     ! H2O
         wmf_tot(49,ilev)=wmf_tot(1,ilev)  ! HDO
         wmf_tot(71,ilev)=wmf_tot(1,ilev)  ! D2O
      end do
c      write(*,*)wmf_tot(1,1),wmf_tot(49,1),wmf_tot(71,1)

c      do k=1,nlev_mav
c        is=specindex(26)+0
c        vpf(is,k)=vpf(is,k)*(1.-2*(1.-t(k)/276)*(t(k)/230.-1.))   ! N2O5
cc
c        is=specindex(31)+0
c        vpf(is,k)=vpf(is,k)*(1.-1.0*(1.-t(k)/296)*(t(k)/220.-1.)) ! CF4
cc
c        is=specindex(33)+0
c        vpf(is,k)=vpf(is,k)*(1.-3*(1.-t(k)/276)*(t(k)/220.-1.))   ! CFC11
cc
c        is=specindex(35)+0
c        vpf(is,k)=vpf(is,k)*(1.-1.8*(1.-t(k)/276)*(t(k)/230.-1.)) ! CCl4
cc
c        is=specindex(42)+0
c        vpf(is,k)=vpf(is,k)*(1.-0.8*(1.-t(k)/276)*(t(k)/230.-1.)) ! HCFC22
cc
c        is=specindex(50)+0
c        vpf(is,k)=vpf(is,k)*(1.-1.5*(1.-t(k)/280)*(t(k)/230.-1.)) ! SF6
cc
c        is=specindex(76)+0
c        vpf(is,k)=vpf(is,k)*(1.-3*(1.-t(k)/280)*(t(k)/220.-1.)) ! CFH2CF3
cc
c      end do

c------------------------------------------------------------------
c  For each species, estimate which layer contains the highest
c  concentration of molecules that are in the ground state.
      do ispeci=1,nspeci
         cgsmax=0.
         lmax=1
c         write(*,*) ' 2: kgas,kiso,ispeci,gasindex(ispeci) = ',
c     &   kgas,kiso,ispeci,gasindex(ispeci)
         do ilev=1,nlev_mav
            vpf(ispeci,ilev)=vpf(ispeci,ilev)*stimem(ilev)*
     &      (trat(ilev)**tdrpf(ispeci))
            cgs=vpf(ispeci,ilev)*d(ilev)*wmf(ispeci,ilev)
            if(abs(cgs).gt.cgsmax) then
               cgsmax=abs(cgs)
               lmax=ilev
            endif
         end do  ! ilev=1,nlev_mav
         maxlev(ispeci)=lmax
c         write(*,*) ispeci,lmax,nlev_mav,cgsmax,wmf(ispeci,lmax)
      end do  ! ispeci=1,nspeci
c---------------------------------------------------------------
c
      tnulst=1.0e-13
c  HITRAN
c      llformat='(i2,i1,f12.6,e10.3,10x,f5.0,f5.4,f10.4,f4.2,f8.6,a24)'
c      llformat='(i2,i1,f12.6,e10.3,10x,2f5.0,f10.4,f4.2,f8.6,a)'
      llformat='(i2,a1,f12.6,e10.3,10x,2f5.0,f10.4,f4.2,f8.6,a)'
      llformat2='(I2,a1,F12.6,E10.3,10x,F10.4,F6.4,F5.3,F9.6,F9.6,F6.4,F
     &5.3,F9.6,F9.6,F5.3,F8.5,F8.5,F8.5,F9.6,F9.6,F9.6,F9.6,F9.6,F9.6,F9
     &.6,F9.6,F9.6,F9.6,F9.6,F9.6,A34)'
c
c Make sure that all linelists are present and correct.
c (Avoids waiting ~15 minutes while abscoj crunches through atm.101
c  only to crash on a subsequent linelist).
c      call substr(linefiles,linfil,mlf,nlf)
c      if(nlf .gt.mlf) stop 'ABSCOF: increase parameter MLF'
      open(19,file=tll_file,status='old')
      atmnv_flag=0
      do lfile=1,mlf  
         read(19,'(34x,a)',end=76) linfil(lfile)
         if(index(linfil(lfile),'atmnv.275').gt.0) atmnv_flag=1
         inquire(file=linfil(lfile),exist=linfil_exist)
         if(linfil_exist.eqv..false.) then
            write(*,*)'Error: '//linfil(lfile)
            stop 'linelist not found'
         end if
         idot=lloc(linfil(lfile),'.')  ! last dot in filename
         read(linfil(lfile)(idot+1:),'(i3)') reclen
         fsib=file_size_in_bytes(lunr_ll,linfil(lfile))
         open(lunr_ll,file=linfil(lfile),access='direct',
     &   form='formatted',status='old',recl=reclen)
         nlines=int(fsib/reclen,kind(reclen))
         if ( nlines*reclen .ne. fsib ) then
            write(*,*)' Linelist size not exactly divisible by ',reclen
            write(*,*)linfil(lfile),fsib
c            stop
         endif
         close(lunr_ll)
      end do !  lfile=1,nlf
      stop 'abscoj: Increase parameter mlf'
76    nlf=lfile-1
      close(19)

c  Do Full line mixing abs calculation if possible
      if(calc_nv.eqv..true.)then
        if(atmnv_flag.le.0) stop
     &  'atmnv.275 linelist not found in telluric_linelists.md5'
        call get_ggg_environment(gggdir,dl)
        lrt=lnbc(gggdir)
c       P-Branch of CH4 2nu3 band
        if(((fzero.gt.5891.06).and.(fzero.lt.5994.15)).or.
     &     ((cen.gt.5891.06).and.(cen.lt.5994.15)).or.
     &     ((fmax.gt.5891.06).and.(fmax.lt.5994.15)).or.
     &     ((fzero.lt.5891.06).and.(fmax.gt.5994.15)))then
c         Create path to the spec data needed
          sfile_path=gggdir(:lrt)//'linelist/ch4_data/'
          sfile='ch4_2nu3_P.dat'
c         Inputs to pass to full LM code
          ispeci=specindex(6)+1
          mm=targmol(ispeci)
c          np=49
          lp=len_trim(sfile_path)
c         Do full line mixing abs calculation
          call CompAbs_qSDVFullLM(sfile_path(:lp),sfile,nlev_mav,t,p,
     &    wmf,d,fzero,grid,ncp,nspeci,ispeci,mm,molewt(ispeci),nu1,nu2,
     &    vvoigt,vac)
        end if
c       Q-Branch of CH4 2nu3 band
        if(((fzero.gt.5998.21).and.(fzero.lt.6004.87)).or.
     &     ((cen.gt.5998.21).and.(cen.lt.6004.87)).or.
     &     ((fmax.gt.5998.21).and.(fmax.lt.6004.87)).or.
     &     ((fzero.lt.5998.21).and.(fmax.gt.6004.87)))then
c         Create path to the spec data needed
          sfile_path=gggdir(:lrt)//'linelist/ch4_data/'
          sfile='ch4_2nu3_Q.dat'
c         Inputs to pass to full LM code
          ispeci=specindex(6)+1
          mm=targmol(ispeci)
c          np=50
          lp=len_trim(sfile_path)
c         Do full line mixing abs calculation
          call CompAbs_qSDVFullLM(sfile_path(:lp),sfile,nlev_mav,t,p,
     &    wmf,d,fzero,grid,ncp,nspeci,ispeci,mm,molewt(ispeci),nu1,nu2,
     &    vvoigt,vac)
        end if
c       R-Branch of CH4 2nu3 band
        if(((fzero.gt.6015.66).and.(fzero.lt.6114.68)).or.
     &     ((cen.gt.6015.66).and.(cen.lt.6114.68)).or.
     &     ((fmax.gt.6015.66).and.(fmax.lt.6114.68)).or.
     &     ((fzero.lt.6015.66).and.(fmax.gt.6114.68)))then
c         Create path to the spec data needed
          sfile_path=gggdir(:lrt)//'linelist/ch4_data/'
          sfile='ch4_2nu3_R.dat'
c         Inputs to pass to full LM code
          ispeci=specindex(6)+1
          mm=targmol(ispeci)
c          np=42
          lp=len_trim(sfile_path)
c         Do full line mixing abs calculation
          call CompAbs_qSDVFullLM(sfile_path(:lp),sfile,nlev_mav,t,p,
     &    wmf,d,fzero,grid,ncp,nspeci,ispeci,mm,molewt(ispeci),nu1,nu2,
     &    vvoigt,vac)
        end if
c       CO2 Q-Branch @ 792 cm-1
        if(((fzero.gt.791.4).and.(fzero.lt.794.3)).or.
     &     ((cen.gt.791.4).and.(cen.lt.794.3)).or.
     &     ((fmax.gt.791.4).and.(fmax.lt.794.3)))then
c         Create path to the spec data needed
          sfile_path=gggdir(:lrt)//'linelist/co2_data/co2_792Q/'
          sfile='co2_792_Q.dat'
c         Inputs to pass to full LM code
          ispeci=specindex(2)+1
          mm=targmol(ispeci)
c          np=27
          lp=len_trim(sfile_path)
c         Do full line mixing abs calculation
          call CompAbs_qSDVFullLM(sfile_path(:lp),sfile,nlev_mav,t,p,
     &    wmf,d,fzero,grid,ncp,nspeci,ispeci,mm,molewt(ispeci),nu1,nu2,
     &    vvoigt,vac)
        end if
      endif !End non-Voigt Full LM
c
c  HERE STARTS THE MAIN LOOP ON THE LINE parameter file.....
      nltnu=0
      do lfile=1,nlf  
         do jspeci=1,nspeci
            fprev(jspeci)=nu1-1.0
            delpl(jspeci)=0.001
         end do
c
c  Set logical variables depending on what kinds of linelists we have
         scia=.false.
         if(index(linfil(lfile),'scia').gt.0) scia=.true.
         fcia=.false.
         if(index(linfil(lfile),'fcia').gt.0) fcia=.true.
         hitflag=.false.
         if(index(linfil(lfile),'hitran').gt.0) hitflag=.true. 

c  Determine number of records and record length in each linelist
         idot=lloc(linfil(lfile),'.')
         read(linfil(lfile)(idot+1:),'(i3)') reclen
         if(reclen.ge.161) then
            lq=94
         else
            lq=reclen-67
         endif
         fsib=file_size_in_bytes(lunr_ll,linfil(lfile))
         open(lunr_ll,file=linfil(lfile),access='direct',
     &   form='formatted',status='old',recl=reclen)
c         nlines=fsib/reclen
         nlines=int(fsib/reclen,kind(reclen))
         llf=lnbc(linfil(lfile))

c  Find record index of first and last lines of relevance.
         kline1=posnall(lunr_ll,nu1,nlines) ! index of the last line with v < NU1
         kline2=posnall(lunr_ll,nu2,nlines) ! index of the last line with v < NU2
c         write(*,*)' Line1,Line2,NLINES =',kline1,kline2,nlines
         write(lunw_rpt,*)' Line1,Line2,NLINES =',kline1,kline2,nlines
         write(lunw_rpt,*)       '    Freq    G  I   S   Stren      E"'
     &   //'   PBHW   Isotopomer         E-wid  Lay  Line#  TG'
c
c         write(*,*) lfile,linfil(lfile),reclen,lq,kline1+1
         do jline=kline1+1,kline2
            if(index(linfil(lfile),'atmnv').gt.0)then
               if(calc_nv.eqv..false.)goto 145 !skip non Voigt 
               read(lunr_ll,llformat2,rec=jline)kgas,ciso,freq,stren,
     &         eprime,abhw,tdpbhw,pshift,tdpshift,sbhw,tdsbhw,spshift,
     &         tdspshift,sdwp,sdshiftp,dn,corr,Yair1,Yair2,Yair3,Yair4,
     &         Yself1,Yself2,Yself3,Yself4,Yh2o1,Yh2o2,Yh2o3,Yh2o4,
     &         quantum(:lq)
                  kiso=ciso2kiso(kgas,ciso)
            else
               read(lunr_ll,llformat,rec=jline) kgas,ciso,freq,stren,
     &         abhw,sbhw,eprime,tdpbhw,pshift,quantum(:lq)
            kiso=ciso2kiso(kgas,ciso)
c            write(*,*)kgas,kiso,freq,stren
c            if(ciso.eq.'B') write(*,*) ciso,kgas,kiso,freq,stren
               if((calc_nv.eqv..true.).and.
     &           (quantum(79:79).eq.'1'))goto 144 !Skip tagged lines 
            endif

c            if(kgas.eq.2 .and. kiso.eq.1) freq=freq*(1.0d0+0.06E-06)

c            write(*,*)jline,kline1,kline2,kgas,kiso,freq
c            pshift=0.0    ! disable pressure shifts
c         if(pshift.ge.0.006) write(*,*)'gas,f,pshift=',kgas,freq,pshift
            kwas=kgas
            if(hitflag) call hitran_to_atmos_gas_numbering(kgas,kiso)
            if(abhw.le.0.0) abhw=0.1
            if(sbhw.le.0.0) sbhw=abhw


c            write(*,*)lfile,jline,kwas,kgas,kiso,freq,stren
c  Check that all linelist gases/species are in isotopomers.dat file
            if(kgas.le.0 .or. kgas.gt.ngas) then
c            write(*,*)'ABSCOJ: Warning: unknown gas in linelist:', kwas
c            write(*,*)'Linelist contains unrecognized gases (NO+, HOBr)'
               cycle
            endif

c  O3 has lots of lines with pshift=0
            if(kgas.eq.3) then
               if(abs(pshift).le.0.0) pshift=-0.001*sngl(freq)/700
            endif

c            if(kgas.eq.7) then
c               tdpbhw=tdpbhw*0.93
c               tdsbhw=tdsbhw*0.93
c               abhw=abhw*1.005
c               sbhw=sbhw*1.005
c            endif

c  Fudge C2H6 Q-branch widths in 2900-3100 cm-1 region
c  to crudely simulate the effects of line mixing
            if(kgas.eq.38) then           ! CH3Cl
               if(abs(freq-825).lt.80) then
                  read(quantum,'(30x,i3)')jlo
                  abhw= 0.061/(1.0+float(jlo-20)/55)
                  sbhw= 1.5625*abhw
c                  write(*,*)freq,stren, eprime, jlo, abhw
               endif
            endif

c  Fudge CH3Cl Q-branch widths in 2900-3100 cm-1 region
c  to crudely simulate the effects of line mixing
            if(kgas.eq.30) then           ! CH3Cl
               if(abs(freq-3000).lt.100) then
c                  write(*,*) linfil(lfile),kgas,kiso,freq,stren,quantum
                  if(reclen.ge.161) then
                     read(quantum,'(30x,2i3,9x,2i3)')jlo,klo,jhi,khi
                     if(jlo.eq.jhi) then     ! Q-branch lines
                        abhw=abhw*0.62/(1.0+0.020*(iabs(klo)-4)) 
                        stren=stren*(1.00+0.032*(iabs(klo)-3))
                     endif
                  endif
               endif
            endif

cc  Fudge strong CO2 widths for sensitivity analysis
c          if(kgas.eq.2 .and. abs(freq-4850).lt.50) then
c             if(index(linfil(lfile),'hitran').gt.0) then  ! HITRAN
c                read(quantum(52:54),'(i3)') jlo
cc                abhw=0.05+0.04/(1.+float(jlo)/20)
c                abhw=0.99*abhw
c                abhw=abhw*(1. - 0.4*float(jlo-15)/(jlo+15))
c             endif
c          endif

c  CO2 Q-branch widths for LM kludge
            if(kgas.eq.2 .and. abs(freq-828.22) .lt. 0.11) then
               abhw=abhw*0.55
            endif

c            if(kgas.eq.2.and.kiso.eq.1.and.abs(freq-980.0).lt.25.) then
c               abhw=abhw*0.99
c            endif

c            if(kgas.eq.2.and.kiso.eq.1.and.abs(freq-3475.0).lt.25.) then
c            if(index(quantum,'  2 1 1 01       0 0 0 01').gt.0) then
c               stren=stren*1.08
c            endif
c            endif

c            if(kgas.eq.2.and.kiso.eq.7.and.abs(freq-3505.0).lt.20.) then
c            if(index(quantum,'  1 0 0 12       0 0 0 01').gt.0) then
c               stren=stren*1.20
c            endif
c            endif

c            if(kgas.eq.2.and.kiso.eq.7.and.abs(freq-3545.0).lt.20.) then
c            if(index(quantum,'  1 0 0 12       0 0 0 01').gt.0) then
c               stren=stren*1.08
c            endif
c            endif


c  NO2 frequency re-calibration for Perrin 15NO2 linelist.
            if(kgas.eq.10 .and. kiso.eq.2) then
               freq=freq+0.00097-0.000015*(freq-1582)
            endif

            ispeci=specindex(kgas)+kiso

c            write(*,*)kgas,kiso,ispeci,freq,stren
  
            if(ispeci.gt.specindex(kgas+1)) then
               write(lunw_rpt,'(a,i2,a1,i2,a9,i8,a4,a)')
     &         'Warning: unknown gas/speci: ',kgas,'/',kiso,
     &        'on line',jline,' of ',linfil(lfile)(:lnbc(linfil(lfile)))
               cycle
            endif
c
            lmax=maxlev(ispeci)
c  Estimate equivalent width of line in the precomputed spectral interval 
            alor=abhw*p(lmax)*trat(lmax)**tdpbhw ! HWHM
            sxcgs=sngl(stren*d(lmax)*wmf(ispeci,lmax)*vpf(ispeci,lmax))
            if(scia) sxcgs=sxcgs*d(lmax)*wmf(ispeci,lmax)
            if(fcia) sxcgs=sxcgs*d(lmax)*(1.0-wmf(ispeci,lmax))
            linwid=2*sqrt(alor**2+sngl(grid**2)) ! FWHM (assumes Doppler width = grid)
            eptf=eprime*tfac(lmax)
            if(eptf.gt.85.0) eptf=85.0  ! prevent Inf
            tnu=sxcgs*exp(eptf)/linwid/spi
c  If line is centered outside window adjust its TNU
            del=sngl(dabs(cen-freq)-hw) ! distance beyond edge of window
            if(del.gt.linwid) tnu=tnu*(linwid/del)**2

            if(kiso.eq.0 .and. kgas.ne.2) then
               delpl(ispeci)=freq-fprev(ispeci)
               fprev(ispeci)=freq
            endif

c            write(*,*) 
c            if(kgas.eq.1)write(*,*)freq,kgas,kiso,ispeci,stren,tnu,tnulst
            if(abs(tnu).lt.tnulst) cycle   ! skip the really weak lines
c            write(*,*) freq,kgas,kiso,ispeci,stren,eprime,tnu,tnulst
            nltnu=nltnu+1                  ! count the used lines
c
c            if(kgas.eq.1)write(*,*)'lmax,p,alor=',freq,lmax,p(lmax),alor
c  Write out the details of lines of discernable strength.
            if(abs(tnu).ge.1.e-11) write(lunw_rpt,fmt_rpt)
     &      freq,kgas,kiso,ispeci,stren,eprime,abhw,
     &      speci_id(ispeci),2.0e+08*tnu,lmax,jline,targmol(ispeci)

c  Some linelists use E"=-1 to denote no info. Set to large value
            if(kiso.gt.0 .and. eprime.lt.0.0d0) eprime=999. 

            mm=targmol(ispeci)  ! =1 for 1'st target gas etc.
            if(mm.gt.mmax) mmax=mm
c            write(*,*)freq,kgas,kiso,ispeci,mm,targmol(ispeci)

c            if(kgas.eq.5 ) abhw=0.99*abhw

cc  Fudge the CH4 P-branch lines of the 2nu3 band.
c            if(kgas.eq.6 .and. abs(freq-5900).lt.100 )
c     &       abhw=abhw-0.002-0.01*(eprime/800)**2   ! fudge 20080724

c            if( kgas.eq.2 .and. kiso .eq.1 ) then
c               if(freq.gt.2038 .and. freq.lt.2076.0) then
c                  stren=stren*1.015
c               endif
c            endif
c            if( kgas.eq.2 .and. kiso .eq.2 ) stren=stren*1.02
c            if( kgas.eq.2 ) sbhw=sbhw*0.992
c            if( kgas.eq.2 .and. kiso.eq.1 ) sbhw=sbhw*1.02
c            if( kgas.eq.2 .and. kiso.eq.2 ) abhw=abhw*0.97
c            if( kgas.eq.2 .and. kiso.eq.2 ) sbhw=sbhw*0.97

cc  Next 21 lines pertain to MATMOS sensitivity analysis
c            if( index(linfil(lfile),'atm.161').gt.0 .and. kgas.eq.2
c     &     .and. freq.lt.860.0 ) stren=stren*0.9
c
c            if( index(linfil(lfile),'atm.161').gt.0 .and. kgas.eq.2
c     &     .and. abs(freq-2980).lt.40.0 ) stren=stren*0.45
c
c            if( index(linfil(lfile),'atm.161').gt.0 .and. kgas.eq.2
c     &     .and. kiso.eq.4 .and. abs(freq-2991).lt. 1.4) stren=stren*3.5
c
c            if( index(linfil(lfile),'atm.161').gt.0 .and. kgas.eq.2
c     &     .and. kiso.eq.5 .and. abs(freq-1865).lt.30.0 ) stren=stren*5
c
c            if( index(linfil(lfile),'atm.161').gt.0 .and. kgas.eq.2
c     &     .and. kiso.eq.6 ) stren=stren*0.20
c
c            if( index(linfil(lfile),'atm.161').gt.0 .and. kgas.eq.2
c     &     .and. kiso.ge.10) stren=stren*0.020
c
c            if( index(linfil(lfile),'atm.161').gt.0 .and. kgas.eq.2
c     &     .and. kiso.eq.11) stren=stren*12

c  Check that pseudolines (kiso=0) are evenly spaced.
            if(kiso.eq.0 .and. kgas.ne.2) then
               if(delpl(ispeci).le.0) then
                  write(*,*)kgas,kiso,freq,stren,fprev(ispeci)
                  write(*,*)'Fatal Error.'
                  write(*,*)' Multiple pseudo-lines at same frequency?'
                  write(*,*)
     &           'Linelist contains species not in isotopologs.dat?'
                  stop
               endif
            endif
c
c  Loop over atmospheric levels
            do i=1,nlev_mav
               if(p(i).gt.4.0) then
                  write(lunw_rpt,*)' Warning: P > 4 atm'
               else
C              Calculate pressure shift with extra parameters
               if(index(linfil(lfile),'atmnv').gt.0)then
C                 SHIFT = p(i)*((pshift+tdpshift*(t(i)-296))*(1-
C     &                   wmf_tot(kgas,i)) + (spshift+
C     &                   tdspshift*(t(i)-296))*wmf_tot(kgas,i))
                   SHIFT = p(i)*((pshift)*(1-wmf_tot(kgas,i)) + 
     &                    (spshift)*wmf_tot(kgas,i))
               endif   
                  vcent=(freq-fzero+dble(pshift*p(i)))/grid
                  if(kiso.eq.0 .and. kgas.ne.2) then
                     dopphwem=sngl(delpl(ispeci)) !set DOPPHWEM to line spacing
                  else
                     dopphwem=4.3014e-7*sngl(freq)*
     &               sqrt(t(i)/molewt(ispeci))
                  endif

C              Calculate the first order Line-mixing parameter
               if(index(linfil(lfile),'atmnv').gt.0)then
                  YT=p(i)*((1-wmf(ispeci,i))*(Yair1*(trat(i)**3)
     &             +Yair2*(trat(i)**2)+Yair3*(trat(i))+Yair4)
     &             +wmf(ispeci,i)*(YSelf1*(trat(i)**3)
     &             +YSelf2*(trat(i)**2)+YSelf3*trat(i)+Yself4)
     &             +wmf(specindex(1)+1,i)*(Yh2o1*(trat(i)**3)
     &             +Yh2o2*(trat(i)**2)+Yh2o3*trat(i)+Yh2o4))
               endif

c  pbhw = abhw*(1-wmf) + sbhw*wmf
                  if(kgas.eq.7) then
c                    write(*,*) 'Heree',i,wmf(ispeci,i),wmf_tot(kgas,i)
                     pbhw=abhw+(wmf_tot(kgas,i)-0.21)*(sbhw-abhw)/0.79
                  elseif(kgas.eq.41) then
                     pbhw=abhw+(wmf_tot(kgas,i)-0.79)*(sbhw-abhw)/0.21
                  else
                     pbhw=abhw+wmf_tot(kgas,i)*(sbhw-abhw)
                  endif
c                  if(kgas.eq.1) write(62,'(i3,12f12.6)')i,p(i),t(i),
c     &   freq,stren,abhw,sbhw,wmf(specindex(1)+1,i),wmf_tot(kgas,i),pbhw

c  Add extra H2O-broadening, assumed to be greater than the air-broadening.
c                  if(kgas.ne.1 .and. kgas.ne.49 .and. kgas.ne.71) 
c     &            pbhw=pbhw+0.80*abhw*abs(wmf(specindex(1)+1,i))
                  abhw_h2o = ewvb(ispeci)*abhw*wmf_tot(1,i)
                  pbhw=pbhw+abhw_h2o

c  Add extra CO2-broadening, assumed to be 1.15 times the air-broadening.
c  This is only important for Venus and Mars.
                  if(kgas.ne.2) pbhw=pbhw+
     &            .15*abhw*abs(wmf(specindex(2)+1,i))

c Now compute the absorption cross-section of this particular line.
               ! Pressure broadening 
               if(index(linfil(lfile),'atmnv').gt.0)then
                  if((kgas.eq.7).and.(kiso.eq.1))then
                    !only use abhw since sbhw not measured
                     HWT=p(i)*(abhw+abhw_h2o)*trat(i)**tdpbhw
                  else
                    HWT=p(i)*(((abhw+abhw_h2o)*trat(i)**tdpbhw)*
     &              (1-wmf_tot(kgas,i))
     &              +(sbhw*trat(i)**tdsbhw)*wmf_tot(kgas,i))
                  endif
                  y = sngl(HWT)/dopphwem
               else
                  y=pbhw*p(i)*trat(i)**tdpbhw/dopphwem !in DW for Voigt
c                  if(kgas.eq.1) write(*,*)'y= ',y,pbhw,abhw,sbhw,
c     &            wmf_tot(kgas,i)
               endif  
                  sxcgs=sngl(stren*d(i)*wmf(ispeci,i)*vpf(ispeci,i))
                  if(scia) sxcgs=sxcgs*d(i)*wmf(ispeci,i)
                  if(fcia) sxcgs=sxcgs*d(i)*(1-wmf(ispeci,i))
                  eptf=eprime*tfac(i)
                  if(eptf.gt.85.0) eptf=85.0  ! prevent overflow of EXP
                  sxcgsorpidw=sxcgs*exp(eptf)/dopphwem/srpi
c                  write(*,*)'sg=',i,sxcgs,eptf,exp(eptf),dopphwem
                  godw=grid/dopphwem   ! primative grid point spacing (Doppler widths)
                  flinwid=dble((3+amin1(y/frac,
     &            sqrt(abs(y*sxcgsorpidw/srpi/tnulst))))/godw)
c                  if(kgas.eq.2) flinwid=dmin1(flinwid,21950.d0) ! Sub-Lorentzian CO2
                  kv1=max0(1,1+int(vcent-flinwid)) ! Start index (primative grid points)
                  kv2=min0(int(vcent+flinwid),ncp) ! End index (primative grid points)
                  nv=kv2-kv1+1    ! Number of primative grid points encompassed by line
c                  if(kgas.eq.1 .and. nv.lt.1) write(*,*)
c     &            'freq,kv1,kv2,nv,flinwid= ',freq,kv1,kv2,nv,flinwid
                  if(nv.gt.0) then
                     x1=godw*(dble(kv1)-vcent)   ! starting point (doppler widths)
                     if(kgas.eq.2) then
                        sublor= 9.0/dopphwem    ! CO2  9.5 cm-1
                     elseif(kgas.eq.1) then
c  For H2O the shape turns sub-Lorentzian very far from line center
                        sublor=5000.0/dopphwem  ! H2O  5000 cm-1
c                        sublor=50.0/dopphwem  ! H2O  50 cm-1
                     else
                        sublor=100.0/dopphwem   ! other 100 cm-1
                     endif
c                     sublor=150.0/dopphwem  ! for ooofti LM test
c  When x=sublor, godw.n = 100/dopphwem
c  When x=sublor, n.grid  = 100 cm-1
                  if(index(linfil(lfile),'atmnv').gt.0)then
                     SDHWT = sdwp*HWT !SD width parameter
                     SDSHIFT = sdshiftp*SHIFT !SD p-shift parameter
                     sg=(dble(kv1)*grid)+fzero !starting grip point cm-1
                     call qSDV(freq,dopphwem,HWT,SDHWT,SHIFT,
     &                   SDSHIFT,sg,YT,nv,grid,vvoigt(kv1)) 
                  else
                     call gct2_humlik(nv,x1,godw,y,sublor,vvoigt(kv1))
                  endif
                     do jv=kv1,kv2  ! Innermost loop over line profile
                        vac(jv,i,mm)=vac(jv,i,mm)-sxcgsorpidw*vvoigt(jv)
                     end do
                  endif   ! nv.gt.0
               endif   ! p(i).gt.1.2
            end do     ! i=1,nlev_mav  Loop over levels
144       continue   ! Skip calulation of abs for a line
         end do       !  do jline=kline1+1,kline2
145      continue     ! Skip a line list
         close(lunr_ll)
      end do         ! end of loop on line files
      write(lunw_rpt,*)nltnu,' lines (all linelists) with S > TNULST'
      close(lunw_rpt)

c  Replace internally-calculated absorption coefficients by externally calculated ones.
      lr=lnbc(runlab)
      lmfile=''
      do jtarg=1,mmax
         if(targmol(specindex(2)+1).eq.jtarg) then
            lmfile=runlab(:12)//'_co2'
            exit
         endif
         if(targmol(specindex(6)+1).eq.jtarg) then
            lmfile=runlab(:12)//'_ch4'
            exit 
         endif
         if(targmol(specindex(7)+1).eq.jtarg)then
            lmfile=runlab(:12)//'_o2'
            exit
         endif
      end do
c      write(*,'(a)') lmfile
c      write(*,*) lnbc(lmfile)
      write(lmfile(lnbc(lmfile)+1:),'(a1,i5.5,a1,i5.5,a7)')
     & '_',nint(fzero+grid),'_',nint(fzero+ncp*grid),'_lm.bin'
c      write(*,*)'lmfile=',lmfile
      inquire(file=lmfile,exist=lmfile_exist)
      if((lmfile_exist.eqv..false.).and.(calc_lm.eqv..true.)) then
         write(*,*)'LM requested but .bin file not found: '//lmfile
         write(*,*)'Please run compute_lm_absco to generate lmfiles.'
         stop
      elseif((lmfile_exist.eqv..true.).and.(calc_lm.eqv..false.)) then
         write(*,*)'Warning: Not using LM file ',lmfile(:lnbc(lmfile))
      endif
      if(lmfile_exist.and.calc_lm) then
c         write(*,*) ' Calling FSIB:'//lmfile
         fsib=file_size_in_bytes(lunr_lm_bin,lmfile)
         write(*,*)' lmfile found: ',lmfile(:lnbc(lmfile))
         if(fsib.ne.12*ncp*nlev_mav) then
            write(*,*) lmfile(:lnbc(lmfile))//' has size: fsib = ',fsib
            write(*,*) 'ncp,nlev_mav = ',ncp,nlev_mav
            write(*,*) 'Expected size = 4 bytes * 3 columns * ncp * nlev
     &_mav = ',4*3*ncp*nlev_mav
            stop 'fsib .ne. 4*3*ncp*nlev_mav '
         endif
         open(lunr_lm_bin,file=lmfile,access='direct',status='old',
     &   recl=12,form='unformatted')
         irec=0
         do ilev=1,nlev_mav
c           write(*,*)' ncp, ilev, nlev_mav =',ncp,ilev,nlev_mav
            do j=1,ncp
               irec=irec+1
               read(lunr_lm_bin,rec=irec)AbsV,AbsY,AbsW
c               vac(j,ilev,jtarg)=-AbsV-AbsY-AbsW             ! CH4: Pine Voigt + LM 
c               vac(j,ilev,jtarg)=vac(j,ilev,jtarg)-AbsY-AbsW ! O2: Hartnam LM & CIA 
c               vac(j,ilev,jtarg)=vac(j,ilev,jtarg)-AbsY      ! O2: Hartnam LM 
               vac(j,ilev,jtarg)=vac(j,ilev,jtarg)-AbsY       ! first order LM only
            end do
         end do
         close(lunr_lm_bin)
c      else
c         write(*,*)' lmfile missing: ',lmfile
      endif   !  (lmfile_exist)

cc  Code to write out monochromatic volume absorption coefficients (if necessary)
c      do jtg=1,9   ! Loop over target gases
c      do ispeci=1,nspeci
c         if(targmol(ispeci).eq.jtg) go to 177
c      end do
c      return
c177   vacfile='vac_tgxx.out'
c      write(vacfile(7:8),'(i2.2)') jtg
c      open(lunw_vac,file=vacfile,status='unknown')
c      write(lunw_vac,*)2,nlev_mav
c      write(lunw_vac,'(a6,80i12)')' Freq ',(ilev-1,ilev=2,nlev_mav)
c      write(*,*) jtg,ispeci,nlev_mav,ncp
c      do jcp=1,ncp   !  Loop over frequencies
c         write(lunw_vac,'(f14.6,8e12.5)')fzero+jcp*grid,
c     &   (-vac(jcp,ilev,jtg)/wmf(ispeci,ilev)/d(ilev),ilev=2,nlev_mav)
c      end do
c      close(lunw_vac)
c      end do   !  Loop over target gases

      return
      end
