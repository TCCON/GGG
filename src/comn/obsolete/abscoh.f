      subroutine abscoh(nlev,t,p,d,nspeci,targmol,vmr,vpf,
     & llsize,linefiles,parfile,fzero,grid,ncp,vac,vv)
c
c  Performs line-by-line calculations of the volume absorption coefficients,
c  which are stored in a 3-D array (frequency, altitude, target gas).
c  First it reads "molparam.dat" and pre-computes the partition functions.
c  Then it reads the linelists and does the line-by-line calculations.
c
c  INPUT Parameters:
c    nlev             I*4   Number of atmospheric levels.
c    T(nlev)          R*4   Temperature Profile (Kelvin)
c    P(nlev)          R*4   Pressure Profile (Atmospheres)
c    D(nlev)          R*4   Number Density Profile (molec.cm-2.km-1)
c    vmr(nspeci,nlev) R*4   Array of vmr profiles for each gas (initially)
c    vpf(nspeci,nlev) R*4   Workspace Array (to store partition functions)
c    nspeci           I*4   Number of different species/isotopomers in PARFILE
c    targmol(nspeci)  I*4   Group to which absorption of each specie belongs.
c    linefiles        C**   Spectral linelists to be read.
c    parfile          C**   File of molecular parameters (e.g. molparam.dat).
c    fzero            R*8   Frequency immediately prior to the first grid point.
c    grid             R*8   Primitive grid point spacing (cm-1).
c    ncp              I*4   Number of primitive spectral grid points.
c    vv(ncp)          R*4   Work Array used to store Voigt lineshape
c
c  OUTPUT Parameter:
c    vac(ncp,nlev,*)  R*4   3-D matrix of volume absorption coefficients
c
c
c  VAC contains the Volume Absorption Coefficient spectrum for each level
c  for each target gas, and so its dot product with the vector of slant paths
c  yields the optical thickness.
c  Absorption coefficients are stored at frequencies vj=fzero+j.grid  j=1,ncp
c  at nlev different levels, and for different groups of species.
c  Note that in the main program VAC is stored as a 1-D array and so we don't
c  have to worry here about its logical versus physical dimensions.
c
      implicit none
      integer*4
     &   kline1,kline2,jline,nlev,ncp,is,lnbc,llf,
     &   nltnu,   ! total number of lines exceeding TNULST
     &   lunr,    ! parfile & linelists
     &   lunw,    ! "abscof.rpt'
     &   kv1, kv2, ! limits of line profile calculation
     &   nv,jv,   ! number of Voigt point = kv2 - kv1 +1
     &   jspeci,  ! Indexes the different isotopomers
     &   ispeci,  ! Indexes the different isotopomers
     &   nspeci,  ! Number of different species/isotopomers in PARFILE
     &   mspeci,  ! Maximum possible number of species/isotopomers
     &   nvmode,   ! Number of different vibrational modes (ie 3N-6) of each gas
     &   mmode,   ! Maximum possible value of NMODE
     &   mgas,ngas,kgas,    ! Number of different gases found in MOLPAR
c     &   cnss, jss,  ! number of sub-strings
     &   mlev,idot,
     &   mlf,nlf, ! maximum & actual number of linefiles
     & posnall,mm,reclen,nlines,kiso,i,j,k,n,
c     & fstat,statb(13),ierr,
     & lmax,lfile,ilev,istat,fsib,file_size_in_bytes
c
      parameter (mlev=151,mgas=64,mlf=10,mmode=30,lunr=21,lunw=20,
     & mspeci=230)
c MISO increased to 230 for Justus' pseudo-gases used for profiling. GCT 19-8-98
      integer*4
     &   targmol(nspeci),
     &   specindex(mgas+1), ! index of species
     &   maxlev(mspeci),    ! level at which each species has its max density
     &   dgen(mmode),       ! degeneracy of vibrational modes
     &   molewt(mspeci)     ! Molecular Weight
c
      real*8
     &   fzero,        ! frequency immediately prior to first grid point
     &   fmax,         ! fzero+ncp*grid
     &   stren,        ! line strength
     &   hw,           ! half-width (cm-1) of precomputed interval
     &   nuoff,        ! search linelist NUOFF cm-1 beyond ends of VAC array 
     &   nu1,nu2,      ! start & end frequencies (cm-1) for linelist search
     &   grid,         ! primitive point spacing (cm-1)
     &   freq,         ! line center frequency (cm-1)
     &   fprev(mspeci),! line center frequency of previous line of JGAS
     &   vcent,        ! line center position in primitive grid points
     &   cen,          ! central frequency of window (cm-1)
     &   x1,           ! start frequency  (for Voigt evaluation)
     &   godw,         ! GRID / Doppler-width
     &   flinwid       ! half-width of computed lineshape
c
      real*4
     &   vv(ncp),       ! Array of Voigt profile
     &   y,            ! Lorentz-width/Doppler-width (used by VOIGT(x,y))
     &   tnu,tnulst,   ! line-center absorbtance & threshold
     &   linwid,       ! approximate line width (cm-1)
     &   del,          ! line center distance beyond edge of VAC interval
     &   pi,           ! 3.14159265...
     &   srpi,         ! SQRT(3.14159265...)
     &   fia,delta,epsilon,    ! Fractional Isotopomeric Abundance
     &   eprime,       ! ground state energy (cm-1)
     &   abhw,         ! air-broadened half-width
     &   sbhw,         ! self-broadened half width
     &   pbhw,         ! pressure-broadened half width
     &   cgs,cgsmax,   ! concentration of ground-state molecules
     &   pshift,       ! line-center frequency pressure shift (cm-1/atm)
     &   trat(mlev),   ! trat(k)=296./t(k)
     &   tfac(mlev),   ! tfac(k)=conx*(trat(k)-1.0)
     &   stimem(mlev), ! stimem(k)=(1.-exp(sf*trat(k)))/(1.-exp(sf))
     &   tdrpf(mspeci),! T-Dependence of Rotational Partition Function
     &   tdpbhw,       ! T-Dependence of Pressure-Broadened Half-Width
     &   vibfrq(mmode),! Array of vibrational frequencies
     &   conx,         ! conx=-1.43881d0/296  ( = hc/kT;  T=296)
     &   dopphwem,     ! Exact HWEM Doppler width
     &   alor,         ! actual Lorentz width (cm-1). Pressure * PBHW
     &   frac,
     &   vibpf,        ! Vibrational Partition Function.
     &   dx,xx,ss,
     &   eta,          ! Sub- super-Lorentzian coefficient (2.0=Lorentzian)
     &   sf,           ! scratch variable used in computation of STIMEM
     &   sxcgs,        ! scratch variable used to save STREN*CGS
     &   sxcgsopidw,   ! scratch variable multiplying VOIGT in innermost loop
     &   sxcgsorpidw,  ! scratch variable multiplying VOIGT in innermost loop
     &   v296,         ! 296 K vibrational partition function
     &   vmr(nspeci,nlev),vpf(nspeci,nlev),
     &   vac(ncp,nlev,0:*),t(*),p(*),d(*)
      parameter (pi=3.1415927,conx=-(1.43881/296))
c
      character
     &   llformat*49,      ! FORMAT statement for linelist
     &   quantum*34,       ! transition quantum numbers
     &   gasname*8,        ! Gas name
     &   linfil(mlf)*100,  ! Name of currently open linelist
c     &   llname*12,       ! Name of linelist in llsize.dat
     &   qq*2,             ! 
     &   parfile*(*),     ! Name of MOLPAR file (usually "molparam.dat")
c     &   pars(nss)*9,     ! parameters to fit
     &   speci_id(mspeci)*24, ! chemical formulae of the species/isotopomers
     &   llsize*(*),      ! String of names of llsize.dat file
     &   linefiles*(*)    ! String of names of linelists
c
      logical hitran, scia,fcia
c====================================================================
      srpi=sqrt(pi)
c      write(*,*)'abscog '
      if(nlev.gt.mlev) stop 'ABOCOF: Increase parameter MLAY'
      if(nspeci.gt.mspeci) stop 'ABOCOF: Increase parameter MSPECI'
      fmax=fzero+ncp*grid
      open(lunw,file='abscof.rpt',status='unknown')
      write(lunw,43)' Pre-computing VACs from ',
     &fzero+grid,'  to ',fmax,' cm-1 :',ncp,
     &' grid points at a primitive point spacing of ',grid, ' cm-1'
 43   format(a,f10.4,a,f10.4,a,/,i6,a,f10.8,a)
c
c   SET PARAMETERS FOR TRANSMISSION CALCULATION
      hw=grid*(ncp-1)/2
      cen=fzero+grid+hw
      frac=0.0001   !  get 99.99% of line absorption
      nuoff=10.D0*dsqrt(2+dble(grid)*ncp*dble(p(2)))   ! see 1992 Toon MEMO
      
      nu1=fzero+grid-nuoff
      nu2=fmax+nuoff
      write(lunw,*)'NU1 =',nu1,'    NU2 =',nu2,'    NUOFF = ',nuoff 
c
c  PRE-CALCULATE TRAT, TFAC, & STIMULATED EMISSION TERM FOR EACH LAYER
c  Assume that the central window wavenumber CEN is close to the line.
      sf=conx*sngl(cen)
      do k=1,nlev
        trat(k)=296./t(k)
        tfac(k)=conx*(trat(k)-1.0)
        stimem(k)=(1.-exp(sf*trat(k)))/(1.-exp(sf))
      end do
c-------------------------------------------------------------------
c  CALCULATE THE VIBRATION PARTITION FUNCTIONSPECIS
      open(unit=lunr,file=parfile,status='old')
c      open(unit=lunr,file=parfile(:lnbc(parfile)),status='old')
      do jspeci=1,nspeci
         call read_isotop(lunr,kgas,kiso,gasname,speci_id(jspeci),
     &   fia,delta,epsilon,molewt(jspeci),tdrpf(jspeci),
     &   vibfrq,dgen,nvmode,mmode,istat)
         if(istat.ne.0) stop 'READ_ISOTOP: ISTAT.NE.0'
         if(kgas.lt.0) stop 'KGAS<0'
         if(kiso.lt.0) stop 'KISO<0'
         specindex(kgas)=jspeci-kiso

c         call lowercase(gasname)
c         targmol(jspeci)=0
c         do jss=1,nss
c            if(gasname//' '.eq.pars(jss)) targmol(jspeci)=jss
c            if(char(kiso+48)//gasname.eq.pars(jss)) targmol(jspeci)=jss
c            if('t'//gasname.eq.pars(jss)) targmol(jspeci)=jss
c         end do
c
c  CALCULATE VPF AT 296K
c        v296=1.0
c        do j=1,nvmode
c          vibfrq(j)=conx*vibfrq(j)
c          v296=v296*(1.-exp(vibfrq(j)))**dgen(j)
c        end do
        v296=vibpf(296.0,vibfrq,dgen,nvmode)  ! VPF at 296 K
c
c  CALCULATE VPFs AT nlev OTHER TEMPERATURES AND DIVIDE BY v296
        do k=1,nlev
c          vpf(jspeci,k)=1.0/v296
c          do j=1,nvmode
c            vpf(jspeci,k)=vpf(jspeci,k)*
c     &      (1.-exp(vibfrq(j)*trat(k)))**dgen(j)
c          end do
          vpf(jspeci,k)=vibpf(t(k),vibfrq,dgen,nvmode)/v296  ! VPF(T)/VPF(296)
        end do
      end do  ! jspeci=1,nspeci
      close(lunr)
      ngas=kgas
      specindex(ngas+1)=jspeci
      if(ngas.gt.mgas) STOP ' ABSCOF: Increase parameter MGAS'


c  Fudge T-dependence of N2O5, CFC-11, CFC-12,  & HCFC-22.
      do k=1,nlev
c        vpf(26,k)=vpf(26,k)*(1.-2*(1.-t(k)/276)*(t(k)/230.-1.))
        is=specindex(26)+0
        vpf(is,k)=vpf(is,k)*(1.-2*(1.-t(k)/276)*(t(k)/230.-1.))     ! N2O5
c
c        vpf(31,k)=vpf(31,k)*(1.-1.0*(1.-t(k)/296)*(t(k)/220.-1.))
        is=specindex(31)+0
        vpf(is,k)=vpf(is,k)*(1.-1.0*(1.-t(k)/296)*(t(k)/220.-1.))   ! CF4
c
c        vpf(33,k)=vpf(33,k)*(1.-3*(1.-t(k)/276)*(t(k)/220.-1.))
        is=specindex(33)+0
        vpf(is,k)=vpf(is,k)*(1.-3*(1.-t(k)/276)*(t(k)/220.-1.))     ! CFC11
c
c        vpf(35,k)=vpf(35,k)*(1.-1.8*(1.-t(k)/276)*(t(k)/230.-1.))
        is=specindex(35)+0
        vpf(is,k)=vpf(is,k)*(1.-1.8*(1.-t(k)/276)*(t(k)/230.-1.))   ! CCl4
c
c        vpf(42,k)=vpf(42,k)*(1.-0.8*(1.-t(k)/276)*(t(k)/230.-1.))
        is=specindex(42)+0
        vpf(is,k)=vpf(is,k)*(1.-0.8*(1.-t(k)/276)*(t(k)/230.-1.))   ! HCFC22
c
c        vpf(50,k)=vpf(50,k)*(1.-1.5*(1.-t(k)/280)*(t(k)/230.-1.))
        is=specindex(50)+0
        vpf(is,k)=vpf(is,k)*(1.-1.5*(1.-t(k)/280)*(t(k)/230.-1.)) ! SF6
c
      end do
c------------------------------------------------------------------
c  For each species, determine which layer contains the highest
c  concentration of molecules that are in the ground state.
      do ispeci=1,nspeci
         cgsmax=0.
         lmax=1
         do ilev=1,nlev
            vpf(ispeci,ilev)=vpf(ispeci,ilev)*stimem(ilev)*
     &      (trat(ilev)**tdrpf(ispeci))
            cgs=vpf(ispeci,ilev)*d(ilev)*vmr(ispeci,ilev)
            if(abs(cgs).gt.cgsmax) then
               cgsmax=abs(cgs)
               lmax=ilev
            endif
         end do  ! ilev=1,nlev
         maxlev(ispeci)=lmax
      end do  ! ispeci=1,nspeci
c---------------------------------------------------------------
c
c      tnulst=1.0e-12*sqrt(p(2))
      tnulst=1.0e-12
c  HITRAN
c      llformat='(i2,i1,f12.6,e10.3,10x,f5.0,f5.4,f10.4,f4.2,f8.6,a24)'
      llformat='(i2,i1,f12.6,e10.3,10x,2f5.0,f10.4,f4.2,f8.6,a34)'
c
c Make sure that all linelists are present and correct.
c (Avoids waiting ~15 minutes while abscoh crunches through atm.101
c  before crashing on a subsequent linelist).
      call substr(linefiles,linfil,mlf,nlf)
      if(nlf .gt. mlf) stop 'ABSCOF: increase parameter MLF'
      do lfile=1,nlf  
         idot=index(linfil(lfile),'.')
         read(linfil(lfile)(idot+1:),'(i3)') reclen
         fsib=file_size_in_bytes(lunr,linfil(lfile))
         open(lunr,file=linfil(lfile),access='direct',
     &   form='formatted',status='old',recl=reclen)
         nlines=fsib/reclen
         if ( nlines*reclen .ne. fsib ) then
            write(*,*)' Linelist size not divisible by ',reclen
            write(*,*)linfil(lfile),fsib
            stop
         endif
         close(lunr)
      end do !  lfile=1,nlf
c
c  HERE STARTS A LOOP ON THE LINE parameter file.....
      nltnu=0
      do lfile=1,nlf  
         do jspeci=1,nspeci
            fprev(jspeci)=nu1-1.0
         end do
c
         scia=.false.
         if(index(linfil(lfile),'scia').gt.0) scia=.true.
         fcia=.false.
         if(index(linfil(lfile),'fcia').gt.0) fcia=.true.
         hitran=.false.
         if(index(linfil(lfile),'hitran').gt.0) hitran=.true.

         idot=index(linfil(lfile),'.')
         read(linfil(lfile)(idot+1:),'(i3)') reclen
         fsib=file_size_in_bytes(lunr,linfil(lfile))
         open(lunr,file=linfil(lfile),access='direct',
     &   form='formatted',status='old',recl=reclen)
         nlines=fsib/reclen
         llf=lnbc(linfil(lfile))
c         write(*,'(i2,i7,1x,a)')lfile,nlines,linfil(lfile)(:llf)
          
         kline1=posnall(lunr,nu1,nlines)  ! index of the last line with v < NU1
         kline2=posnall(lunr,nu2,nlines)  ! index of the last line with v < NU2
c         write(*,*)' Line1,Line2,NLINES =',kline1,kline2,nlines
         write(lunw,*)' Line1,Line2,NLINES =',kline1,kline2,nlines
 94      format(a,f9.3,a,2i7)
         write(lunw,*)       '    Freq    G  I  S   Stren      E"'
     &   //'   PBHW   Isotopomer         E-wid  Lay  Line#  TG'
c
         do 7749 jline=kline1+1,kline2
            read(lunr,llformat,rec=jline)
     &      kgas,kiso,freq,stren,abhw,sbhw,eprime,tdpbhw,pshift,quantum
c            write(*,*)jline,kline1,kline2,kgas,kiso,freq
c            pshift=0.0    ! disable pressure shifts
            if(index(linfil(lfile),'hitran').gt.0) then ! Re-map HITRAN gas numbers to the ATMOS 
               call hitran_to_atmos_gas_numbering(kgas,kiso)
c               if(kgas.eq.1 .and. kiso.ge.4) then !    HDO
c                  kgas=49
c                  kiso=kiso-3
c               else if (kgas .eq. 22) then            ! N2
c                  kgas = 41
c               else if (kgas .eq. 23) then            ! HCN
c                  kgas = 28
c               else if (kgas .eq. 24) then            ! CH3Cl
c                  kgas = 30
c               else if (kgas .eq. 25) then            ! H2O2
c                  kgas  = 23
c               else if (kgas .eq. 26) then            ! C2H2
c                  kgas = 40
c               else if (kgas .eq. 27) then            ! C2H6
c                  kgas = 38
c               else if (kgas .eq. 28) then            ! PH3
c                  kgas = 56
c               else if (kgas .eq. 29) then            ! COF2
c                  kgas = 36
c               else if (kgas .eq. 30) then            ! SF6
c                  kgas = 50
c               else if (kgas .eq. 31) then            ! H2S
c                  kgas = 47
c               else if (kgas .eq. 32) then            ! HCOOH
c                  kgas = 46
c               else if (kgas .eq. 33) then            ! HO2
c                  kgas = 22
c               else if (kgas .eq. 34) then            ! Atomic O
c                  kgas = 00
c               else if (kgas .eq. 35) then            ! ClNO3
c                  kgas = 27
c               else if (kgas .eq. 36) then            ! NO+
c                  kgas = 00
c               else if (kgas .eq. 37) then            ! HOBr
c                  kgas = 00
c               else if (kgas .eq. 38) then            ! C2H4
c                  kgas = 39
c               endif
            endif

cc  This code allowd the CO2 continuum absorption to the fitted separately
cc  from the discrete CO2 lines in the 2400 cm-1 region
c        if(kgas.eq.2) then
c           if(abs(freq-2350).lt.35) then
c              kiso=8
c           endif
c        endif
c
c  The following code block was commented out 2006-03-06
c  O2 super-lorentzian far-wing lineshape
c           eta=2.0
c           if(kgas.eq.7 .and. kiso.ne.0) then
c              if(abs(freq-7885).le.250.)then 
c                 eta=1.986
c                 eta=1.972
c              elseif( abs(freq-13082).le.250.)then 
c                 eta=1.977
c              endif
c           endif
c
cc  The following code renormalizes Aaron Goldman's O2 line strengths
cc  so that they are the same as Lafferty's in the high-J limit.
c          if(kgas.eq.7 .and. kiso.ne.0 .and. abs(freq-7885).le.250.)then 
c             read(quantum,'(15x,a1,i2,a1,i2)') qq(1:1),n,qq(2:2),j
c             if(qq.eq.'OP') stren=stren*(0.93-0.0005*n)
c             if(qq.eq.'PP') stren=stren*(0.97)
c             if(qq.eq.'PQ') stren=stren*(0.97-0.0005*n)
c             if(qq.eq.'QP') stren=stren*1.05
c             if(qq.eq.'QQ') stren=stren*1.00
c             if(qq.eq.'QR') stren=stren*0.97
c             if(qq.eq.'RQ') stren=stren*(1.05+0.0010*n)
c             if(qq.eq.'RR') stren=stren*1.05
c             if(qq.eq.'SR') stren=stren*(1.05+0.0025*n)
c             stren=stren*0.89  ! fudge strengths to give right answer
cc             abhw=abhw*1.04    ! fudge widths to remove "smile"
cc             tdpbhw=tdpbhw*1.05 ! fudge width to remove T-dependence
c          endif
c
cc  Replace O2 widths with calculation based on fits to laboratory O2
cc  spectra of Brown et al, Yang et al, and Ritter & Wilkinson (Mar 2002).
c          if(kgas.eq.7 .and. kiso.ne.0 .and. abs(freq-7885).le.250.) then 
c             read(quantum,'(15x,a1,i2,a1,i2)') dn,n,dj,j
c             if(dn.eq.'O') delta_n=-2
c             if(dn.eq.'P') delta_n=-1
c             if(dn.eq.'Q') delta_n=-0
c             if(dn.eq.'R') delta_n=+1
c             if(dn.eq.'S') delta_n=+2
c             nprime=n+delta_n
c             zz=float(nprime-3*delta_n)
c             pshift=-(8.0E-07*freq*(zz+5)/(zz+15))
c             xx=float(nprime)
c             sbhw=0.02204+0.03749/
c     &       (1+0.05428*xx-0.119E-02*xx**2+0.2073E-05*xx**4)
c             abhw=sbhw*1.012/sqrt(1+((xx-2.5)/45)**2)
c             tdpbhw=0.76/sqrt(1+((xx-5)/55)**2)
c          endif
c          endif
c
c
cc  Reduce CH4 line widths in 6000 cm-1 region.
cc  This code must be removed when the CH4 linelist with Rebecca's
cc  new widths is used.
c          if(kgas.eq.6) then
c             if(abs(freq-5930).le.220.) abhw=abhw*0.97
c          endif

          if(kgas.eq.10) then  ! Set NO2 widths to Dana (JQSRT 57, 445-457, 97)
            tdpbhw=.75
c            read(quantum,'(15x,i2)')nprime2
c            abhw=0.001*( 84.1 - .753*nprime2 + .0059*nprime2**2)
          endif
c          if(kgas.eq.23 .and. lfile.eq.2) go to 7749 !  don't use H2O2 lines from atm.h92

    
c  Check that all linelist gases/species are in isotopomers.dat file
          if(kgas.le.0 .or. kgas.gt.ngas) then
             Write(*,*)'ABSCOG: Warning: unknown gas in linelist:',
     &       kgas
             go to 7749
          endif

          ispeci=specindex(kgas)+kiso

          if(ispeci.gt.specindex(kgas+1)) then
             Write(*,*)'ABSCOG: Warning: unknown speci:',kgas,kiso,
     &       ' on line ',jline,' of ',linfil(lfile)
             go to 7749
          endif
c
          lmax=maxlev(ispeci)
c  Estimate equivalent width of line in the precomputed spectral interval 
          alor=abhw*p(lmax)*trat(lmax)**tdpbhw ! HWHM
          sxcgs=stren*d(lmax)*vmr(ispeci,lmax)*vpf(ispeci,lmax)
          if(scia) sxcgs=sxcgs*d(lmax)*vmr(ispeci,lmax)
          if(fcia) sxcgs=sxcgs*d(lmax)*(1.0-vmr(ispeci,lmax))
c          linwid=2*sqrt(alor**2+dopphwem**2) ! FWHM (assumes Doppler width = grid)
          linwid=2*sqrt(alor**2+grid**2) ! FWHM (assumes Doppler width = grid)
          tnu=sxcgs*exp(eprime*tfac(lmax))/linwid/pi
c  If line is centered outside window adjust its TNU
          del=sngl(dabs(cen-freq)-hw) ! distance beyond edge of window
          if(del.gt.linwid) tnu=tnu*(linwid/del)**2
          if(abs(tnu).lt.tnulst) go to 7748
          nltnu=nltnu+1
c
          if(abs(tnu).ge.2.e-10) then
            write(lunw,7764)freq,kgas,kiso,ispeci,stren,eprime,abhw,
     &      speci_id(ispeci),2000*1.0e+05*tnu,lmax,jline,targmol(ispeci)
          endif

 7764   format(f10.3,i3,2i3,1pe10.2,0pf8.1,1x,f4.3,2x,a18,f7.1,i4,i8,i4)
c
          if(kiso.gt.0 .and. eprime.lt.0.0d0) eprime=999.    ! CH4 lines with E"=-1 are common
          mm=targmol(ispeci)  ! =1 for 1'st target gas etc.
c          if(kgas.eq.2 .or. kgas.eq.6 .or. kgas.eq.7) abhw=1.02*abhw   ! fudge
          if(kiso.eq.0) then
             if(freq-fprev(ispeci).le.0) then
c                write(*,*)kgas,kiso,freq,fprev(ispeci)
                stop 'multiple pseudo-lines at same frequency'
             endif
          endif
c
c
          do i=1,nlev
            vcent=(freq-fzero+dble(pshift*p(i)))/grid
c            vcent=(freq-fzero+dble(pshift*p(i)*trat(i)**0.75))/grid  ! 05-12-26
            if(kiso.eq.0) then
              dopphwem=sngl(freq-fprev(ispeci)) !set DOPPHWEM to line spacing
            else
              dopphwem=4.3e-7*sngl(freq)*sqrt(t(i)/molewt(ispeci))
            endif
            if(sbhw.le.0.0) sbhw=abhw
c            pbhw=abhw*(1.-vmr(ispeci,i))) + sbhw*vv
            if(kgas.eq.7) then
               pbhw=abhw+(vmr(ispeci,i)-0.21)*(sbhw-abhw)/(1-0.21)
               pbhw=pbhw*1.00
            elseif(kgas.eq.41) then
               pbhw=abhw+(vmr(ispeci,i)-0.79)*(sbhw-abhw)/(1-0.79)
            else
               pbhw=abhw+vmr(ispeci,i)*(sbhw-abhw)
            endif

c  Add H2O-broadening, assumed to be 1.35 times the air-broadening
            if(kgas.ne.1 .and. kgas.ne.49) pbhw=pbhw+0.35*abhw*vmr(1,i)

c  Fudge widths of first target gas.
c            if(targmol(ispeci).eq.1) pbhw=pbhw*0.95

            y=pbhw*p(i)*trat(i)**tdpbhw/dopphwem
            sxcgs=stren*d(i)*vmr(ispeci,i)*vpf(ispeci,i)
            if(scia) sxcgs=sxcgs*d(i)*vmr(ispeci,i)
            if(fcia) sxcgs=sxcgs*d(i)*(1-vmr(ispeci,i))
c	    if(kiso.eq.0) then  !  O2 and N2 Collision Induced Absorption (CIA)
c                if(kgas.eq.7 .or. kgas.eq.41) sxcgs=sxcgs*p(i)
c            endif
            sxcgsopidw=sxcgs*exp(eprime*tfac(i))/dopphwem/pi
            godw=grid/dopphwem
            flinwid=dble((3+amin1(y/frac,
     &      sqrt(abs(y*sxcgsopidw/tnulst))))/godw)
            if(kgas.eq.2) flinwid=dmin1(flinwid,21950.d0)  ! Sub-Lorentzian CO2
            kv1=max0(1,1+int(vcent-flinwid))
            kv2=min0(int(vcent+flinwid),ncp)
            nv=kv2-kv1+1
            if(nv.gt.0) then
              x1=godw*(dble(kv1)-vcent)  ! start
              call humlik(nv,x1,godw,y,vv(kv1))
              sxcgsorpidw=sxcgsopidw*srpi
              ss=sxcgsorpidw
              do jv=kv1,kv2
c                 if(eta.ne.2.0) then
c                    xx=(jv-vcent)*godw/(y+0.01)
c                    ss=sxcgsorpidw*(1.0+xx**2)**(1.0-eta/2)
c                 endif
c                 dx=-0.0006*(jv-vcent)
c                 vac(jv,i,mm)=vac(jv,i,mm)-
c     &           ss*vv(jv)*(1.+dx/(1+sqrt(1+dx**2)))
                 vac(jv,i,mm)=vac(jv,i,mm)-ss*vv(jv)
              end do
c              call vsma(vv(kv1),1,-sxcgsorpidw,vac(kv1,i,mm),1,
c     &        vac(kv1,i,mm),1,nv)
            endif   ! nv.gt.0
          end do  ! i=1,nlev
 7748     continue
          if(kiso.eq.0) fprev(ispeci)=freq
 7749   continue
        close(lunr)
      end do     ! end of loop on line files
      write(lunw,*)nltnu,' lines (all linelists) with S > TNULST'
      close(lunw)

c      write(25,*) 2,2
c      write(25,*) 'f a'
c      do j=1,ncp
c         write(25,*) fzero+j*grid,-vac(j,1,1)
c      end do
      return
      end

      function vibpf(temp,vibfrq,dgen,nvmode)
c Calculates the vibrational partition function
c the fraction of molecules that are in the vibrational ground state.
      integer*4 imode,nvmode,dgen(nvmode)
      real*4 temp,vibfrq(nvmode),vibpf
      vibpf=1.0
      do imode=1,nvmode
        vibpf=vibpf*(1.-exp(-(1.43881*vibfrq(imode)/temp)))**dgen(imode)
      end do
      return
      end

