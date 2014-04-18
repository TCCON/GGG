      subroutine abscoj(runlab,nlev,t,p,d,nspeci,targmol,vmr,vpf,
     & linefiles,parfile,fzero,grid,ncp,vac,vvoigt)
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
c    D(nlev)          R*4   Number Density Profile (molec.cm-3)
c    nspeci           I*4   Number of different species/isotopomers in PARFILE
c    targmol(nspeci)  I*4   Group to which absorption of each specie belongs.
c    vmr(nspeci,nlev) R*4   Array of vmr profiles for each gas (initially)
c    vpf(nspeci,nlev) R*4   Workspace Array (to store partition functions)
c    linefiles        C**   Spectral linelists to be read.
c    parfile          C**   File of molecular parameters (e.g. molparam.dat).
c    fzero            R*8   Frequency immediately prior to the first grid point.
c    grid             R*8   Primitive grid point spacing (cm-1).
c    ncp              I*4   Number of primitive spectral grid points.
c    vvoigt(ncp)      R*4   Work Array used to store Voigt lineshape
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
      include "../ggg_const_params.f"
      include "../ggg_int_params.f"

      integer*4
     &   kline1,kline2,jline,nlev,ncp,is,lnbc,lloc,llf,lr,j,
     &   irec,lq,
     &   jlo,jhi,
     &   klo,khi,
c     &   jtg, jcp,
     &   nltnu,    ! total number of lines exceeding TNULST
     &   lun_ll,   ! parfile & linelists
     &   lun_iso,  ! isotopologs
     &   lunw_rpt,     ! 'abscof.rpt'
     &   lun_vac,  ! 'vac_ilev01.out'
     &   lun_rbin, ! reading binary LM files
     &   kv1, kv2, ! limits of line profile calculation
     &   nv,jv,    ! number of Voigt point = kv2 - kv1 +1
     &   jspeci,   ! Indexes the different isotopomers
     &   ispeci,   ! Indexes the different isotopomers
     &   nspeci,   ! Number of different species/isotopomers in PARFILE
     &   nvmode,   ! Number of different vibrational modes (ie 3N-6) of each gas
     &   ngas,kgas,jgas,  ! Number of different gases found in MOLPAR
     &   idot,
     &   mlf,nlf,  ! maximum & actual number of linefiles
     &   posnall,jtarg,mm,mmax,reclen,nlines,kiso,i,k,
     &   lmax,lfile,ilev,istat

      integer*8 fsib,file_size_in_bytes
c
      parameter (mlf=10,lun_ll=21,lun_iso=22,
     &   lunw_rpt=23, lun_vac=24, lun_rbin=46)
      integer*4
     &   targmol(nspeci),
     &   specindex(mgas+1), ! specindex(kgas)+kiso  is the specie number/identifier
     &   gasindex(nspeci),  ! gasindex(ispeci) is the gas number
     &   maxlev(nspeci),    ! level at which each species has its max density
     &   dgen(mvmode),      ! degeneracy of vibrational modes
     &   molewt(nspeci)     ! Molecular Weight
c
      real*8
     &   fzero,        ! frequency immediately prior to first grid point
     &   fmax,         ! fzero+ncp*grid
     &   fcen,         ! fzero+ncp*grid/2
     &   stren,        ! line strength
     &   hw,           ! half-width (cm-1) of precomputed interval
     &   nuoff,        ! search linelist NUOFF cm-1 beyond ends of VAC array 
     &   nu1,nu2,      ! start & end frequencies (cm-1) for linelist search
     &   grid,         ! primitive point spacing (cm-1)
     &   freq,         ! line center frequency (cm-1)
     &   delpl(nspeci),! frequency spacing of consecutive pseudo-lines
     &   fprev(nspeci),! line center frequency of previous line of JSPECI
     &   vcent,        ! line center position in primitive grid points
     &   cen,          ! central frequency of window (cm-1)
     &   x1,           ! start frequency  (for Voigt evaluation)
     &   godw,         ! GRID / Doppler-width
     &   flinwid       ! half-width of computed lineshape
c
      real*4
     &   yy,
     &   vv,zz,
     &   vvoigt(ncp),  ! Array of Voigt profile
     &   y,            ! Lorentz-width/Doppler-width (used by VOIGT(x,y))
     &   tnu,tnulst,   ! line-center absorbtance & threshold
     &   linwid,       ! approximate line width (cm-1)
     &   del,          ! line center distance beyond edge of VAC interval
     &   srpi,         ! SQRT(3.14159265...)
     &   fia(nspeci),  ! Fractional Isotopolog Abundance
     &   delta,epsilon,! Fractionation at tropopause and its gradient
     &   eprime,       ! ground state energy (cm-1)
     &   eptf,         ! eprime*tfac
     &   abhw,         ! air-broadened half-width
     &   sbhw,         ! self-broadened half width
     &   pbhw,         ! pressure-broadened half width
     &   cgs,cgsmax,   ! concentration of ground-state molecules
     &   pshift,       ! line-center frequency pressure shift (cm-1/atm)
     &   trat(mlev),   ! trat(k)=296./t(k)
     &   tfac(mlev),   ! tfac(k)=conx*(trat(k)-1.0)
     &   stimem(mlev), ! stimem(k)=(1.-exp(sf*trat(k)))/(1.-exp(sf))
     &   tdrpf(nspeci),! T-Dependence of Rotational Partition Function
     &   tdpbhw,       ! T-Dependence of Pressure-Broadened Half-Width
     &   vibfrq(mvmode),! Array of vibrational frequencies
     &   conx,         ! conx=-1.43881d0/296  ( = hc/kT;  T=296)
     &   dopphwem,     ! Exact HWEM Doppler width
     &   sublor,       ! line width after which sub-Lorentzian wings kicks in
     &   alor,         ! actual Lorentz width (cm-1). Pressure * PBHW
     &   frac,
     &   vibpf,        ! Vibrational Partition Function.
     &   sf,           ! scratch variable used in computation of STIMEM
     &   sxcgs,        ! scratch variable used to save STREN*CGS
     &   sxcgsorpidw,  ! scratch variable multiplying VOIGT in innermost loop
     &   v296,         ! 296 K vibrational partition function
     &   vmr_tot(mgas,nlev), ! Total FIA-weighted vmr for each gas
     &   vmr(nspeci,nlev),vpf(nspeci,nlev),vmax,
     &   vac(ncp,nlev,0:*),t(nlev),p(nlev),d(nlev)
      parameter (conx=-(1.43881/296))
c
      character
     & runlab*(*),         ! spectrum name
     & lmfile*120,         ! file containing LM absorption coefficients
c     &   vacfile*48,
c     &   rcsmax*40,
     &   fmt_rpt*80,
     &   llformat*49,      ! FORMAT statement for linelist
     &   quantum*94,       ! transition quantum numbers
     &   gasname*8,        ! Gas name
     &   linfil(mlf)*(mfilepath),  ! Name of currently open linelist
     &   parfile*(*),      ! Name of MOLPAR file (usually "molparam.dat")
     &   speci_id(nspeci)*24, ! chemical formulae of the species/isotopomers
     &   linefiles*(*)     ! String of names of linelists
c
      logical scia,fcia,lmfile_exist
c====================================================================
      fmt_rpt='(f10.3,3i3,1pe10.2,0pf8.1,f5.3,2x,a18,f7.1,i4,i8,i4)'
      mmax=0
      lq=0   ! Avoid compiler warning
      srpi=sqrt(spi)
c      write(*,*)'abscog '
      if(nlev.gt.mlev) then
         write(*,*)'nlev,mlev=',nlev,mlev
         stop 'ABSCO: Increase parameter MLEV'
      endif
c      if(nspeci.gt.mspeci) stop 'ABOCOF: Increase parameter MSPECI'
      fmax=fzero+ncp*grid
      fcen=fzero+ncp*grid/2
      open(lunw_rpt,file='abscof.rpt',status='unknown')
      write(lunw_rpt,43)' Pre-computing VACs from ',
     &fzero+grid,'  to ',fmax,' cm-1 :',ncp,
     &' grid points at a primitive point spacing of ',grid, ' cm-1'
 43   format(a,f12.6,a,f12.6,a,/,i6,a,f12.10,a)
c
c   SET PARAMETERS FOR TRANSMISSION CALCULATION
      hw=grid*(ncp-1)/2
      cen=fzero+grid+hw
      frac=5.0E-06   !  capture 99.9975% of line absorption
      nuoff=10.D0*dsqrt(2+dble(grid)*ncp*dble(vmax(p,1,nlev)))
      if(nuoff.gt.700) nuoff=700.0  ! Venus
c      write(*,*)'pmax=',vmax(p,1,nlev)
c      write(*,*)'nuoff=',nuoff,grid*ncp
      
      nu1=fzero+grid-nuoff
      nu2=fmax+nuoff
      write(lunw_rpt,*)'NU1 =',nu1,'    NU2 =',nu2,'    NUOFF = ',nuoff 
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
c  CALCULATE THE VIBRATION PARTITION FUNCTIONS (VPF)
      open(unit=lun_iso,file=parfile,status='old')
c      open(unit=lun_iso,file=parfile(:lnbc(parfile)),status='old')
      do jspeci=1,nspeci
         call read_isotop(lun_iso,kgas,kiso,gasname,speci_id(jspeci),
     &   fia(jspeci),delta,epsilon,molewt(jspeci),tdrpf(jspeci),
     &   vibfrq,dgen,nvmode,mvmode,istat)
         if(istat.ne.0) stop 'READ_ISOTOP: ISTAT.NE.0'
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
        v296=vibpf(296.0,vibfrq,dgen,nvmode)  ! VPF at 296 K
c
c  CALCULATE VPFs AT nlev OTHER TEMPERATURES AND DIVIDE BY v296
        do k=1,nlev
          vpf(jspeci,k)=vibpf(t(k),vibfrq,dgen,nvmode)/v296  ! VPF(T)/VPF(296)
        end do
      end do  ! jspeci=1,nspeci
      close(lun_iso)
      ngas=kgas
      specindex(ngas+1)=jspeci
      if(ngas.gt.mgas) STOP ' ABSCOF: Increase parameter MGAS'

c  Compute total vmr for each gas (used for self-broadening).
c  This is done by adding the true vmrs (FIA*VMR) for each isotopolog
        do ilev=1,nlev   ! Loop over Levels
           jgas=1
           jspeci=1
           vmr_tot(jgas,ilev)=fia(jspeci)*vmr(jspeci,ilev)
           do jspeci=2,nspeci   ! Loop over Gases
c  Sum over isotopologs of each gas
              if(jgas.eq.gasindex(jspeci)) then
c              write(*,*)jspeci,jgas,gasindex(jspeci)
                 vmr_tot(jgas,ilev)=vmr_tot(jgas,ilev)+
     &           fia(jspeci)*vmr(jspeci,ilev)
              else
                 jgas=jgas+1
                 vmr_tot(jgas,ilev)=fia(jspeci)*vmr(jspeci,ilev)
              end if
c              if(ilev.eq.3) write(*,*) jspeci,jgas,fia(jspeci),
c     &        vmr(jspeci,ilev),vmr_tot(jgas,ilev)
           end do
c  Water vapour is a special case, H2O, HDO & D2O have to be added.
c  And the FIA of D2O are currently those of H2O (format issue)
c  and so must be multiplied by the square of the IUPAC D/H ratio.
           vmr_tot(1,ilev)=vmr_tot(1,ilev)+vmr_tot(49,ilev)+
     &     9.653449e-08*vmr_tot(71,ilev)     ! H2O
           vmr_tot(49,ilev)=vmr_tot(1,ilev)  ! HDO
           vmr_tot(71,ilev)=vmr_tot(1,ilev)  ! D2O
        end do
c        write(*,*)vmr_tot(1,1),vmr_tot(49,1),vmr_tot(71,1)

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
c         write(*,*) ispeci,lmax,nlev,cgsmax,vmr(ispeci,lmax)
      end do  ! ispeci=1,nspeci
c---------------------------------------------------------------
c
      tnulst=1.0e-13
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
         idot=lloc(linfil(lfile),'.')  ! last dot in filename
         read(linfil(lfile)(idot+1:),'(i3)') reclen
         fsib=file_size_in_bytes(lun_ll,linfil(lfile))
         lq=reclen-67
         open(lun_ll,file=linfil(lfile),access='direct',
     &   form='formatted',status='old',recl=reclen)
         nlines=fsib/reclen
         if ( nlines*reclen .ne. fsib ) then
            write(*,*)' Linelist size not divisible by ',reclen
            write(*,*)linfil(lfile),fsib
c            stop
         endif
         close(lun_ll)
      end do !  lfile=1,nlf
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

c  Determine number of records and record length in each linelist
         idot=lloc(linfil(lfile),'.')
         read(linfil(lfile)(idot+1:),'(i3)') reclen
         fsib=file_size_in_bytes(lun_ll,linfil(lfile))
         open(lun_ll,file=linfil(lfile),access='direct',
     &   form='formatted',status='old',recl=reclen)
         nlines=fsib/reclen
         llf=lnbc(linfil(lfile))
          
c  Find record index of first and last lines of relevance.
         kline1=posnall(lun_ll,nu1,nlines) ! index of the last line with v < NU1
         kline2=posnall(lun_ll,nu2,nlines) ! index of the last line with v < NU2
c         write(*,*)' Line1,Line2,NLINES =',kline1,kline2,nlines
         write(lunw_rpt,*)' Line1,Line2,NLINES =',kline1,kline2,nlines
         write(lunw_rpt,*)       '    Freq    G  I  S   Stren      E"'
     &   //'   PBHW   Isotopomer         E-wid  Lay  Line#  TG'
c
         do jline=kline1+1,kline2
            read(lun_ll,llformat,rec=jline) kgas,kiso,freq,stren,
     &      abhw,sbhw,eprime,tdpbhw,pshift,quantum(:lq)
            if(kgas.eq.2 .and. kiso.eq.0) kiso=10  ! HITRAN 2012 kludge
c            write(*,*)jline,kline1,kline2,kgas,kiso,freq
c            pshift=0.0    ! disable pressure shifts
c         if(pshift.ge.0.006) write(*,*)'gas,f,pshift=',kgas,freq,pshift
            if(index(linfil(lfile),'hitran').gt.0) then ! Re-map HITRAN gas numbers to the ATMOS 
               call hitran_to_atmos_gas_numbering(kgas,kiso)
            endif
            if(abhw.le.0.0) abhw=0.1
            if(sbhw.le.0.0) sbhw=abhw

c  Check that all linelist gases/species are in isotopomers.dat file
          if(kgas.le.0 .or. kgas.gt.ngas) then
c          write(*,*)'ABSCOI: Warning: unknown gas in linelist:', kgas
c          write(*,*)'Linelist contains unrecognized gases (NO+, HOBr)'
             cycle
          endif

c  Fudge CH3Cl Q-branch widths in 2900-3100 cm-1 region
c  to crudely simulate the effects of line mixing
c          if(index(linfil(lfile),'ch3cl_apfud').gt.0) then
          if(kgas.eq.30) then           ! CH3Cl
             if(abs(freq-3000).lt.100) then
                if(index(linfil(lfile),'hitran').gt.0) then  ! HITRAN
                   read(quantum,'(30x,2i3,9x,2i3)')jlo,klo,jhi,khi
                else
                   read(quantum,'(8x,2i3,9x,2i3)')jlo,klo,jhi,khi
                endif
                if(jlo.eq.jhi) then     ! Q-branch lines
                   abhw=abhw*0.62/(1.0+0.020*(iabs(klo)-4)) 
                   stren=stren*(1.00+0.032*(iabs(klo)-3))
c                   kiso=3
                endif
             endif
          endif
c          endif

          ispeci=specindex(kgas)+kiso

c          write(*,*)kgas,kiso,ispeci,freq,stren

          if(ispeci.gt.specindex(kgas+1)) then
             write(lunw_rpt,'(a,i2,a1,i2,a9,i8,a4,a)')
     &       ' Warning: unknown gas/speci: ',kgas,'/',kiso,
     &       ' on line',jline,' of ',linfil(lfile)(:lnbc(linfil(lfile)))
             cycle
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

c          write(*,*) 
c          write(*,*) freq,kgas,kiso,ispeci,stren,eprime,tnu,tnulst
          if(abs(tnu).lt.tnulst) cycle   ! skip the really weak lines
c          write(*,*) freq,kgas,kiso,ispeci,stren,eprime,tnu,tnulst
          nltnu=nltnu+1                  ! count the used lines
c
c  Write out the details of lines of discernable strength.
          if(abs(tnu).ge.1.e-11) write(lunw_rpt,fmt_rpt) freq,kgas,kiso,
     &    ispeci,stren,eprime,abhw,
     &       speci_id(ispeci),2.0e+08*tnu,lmax,jline,targmol(ispeci)

c  Some linelists use E"=-1 to denote no info. Set to large value
          if(kiso.gt.0 .and. eprime.lt.0.0d0) eprime=999. 

          mm=targmol(ispeci)  ! =1 for 1'st target gas etc.
          if(mm.gt.mmax) mmax=mm
c          write(*,*)freq,kgas,kiso,ispeci,mm,targmol(ispeci)

c  Fudge the CH4 P-branch lines of the 2nu3 band.
          if(kgas.eq.6 .and. abs(freq-5900).lt.100 )
     &     abhw=abhw-0.002-0.01*(eprime/800)**2   ! fudge 20080724

c  Check that pseudolines (kiso=0) are evenly spaced.
          if(kiso.eq.0 .and. kgas.ne.2) then
            if(delpl(ispeci).le.0) then
               write(*,*)kgas,kiso,freq,stren,fprev(ispeci)
               write(*,*)'Fatal Error.'
               write(*,*)' Multiple pseudo-lines at same frequency?'
               write(*,*)
     &        'Linelist contains species not in isotopologs_local.dat?'
               stop
            endif
          endif
c
c  Loop over atmospheric levels
          do i=1,nlev
             if(p(i).lt.1.1) then
             vcent=(freq-fzero+dble(pshift*p(i)))/grid
c             vcent=(freq-fzero+dble(pshift*p(i)*trat(i)**0.75))/grid  ! 05-12-26
c             if(kiso.eq.0 .and. kgas.ne.32) then
             if(kiso.eq.0 .and. kgas.ne.2) then
               dopphwem=sngl(delpl(ispeci)) !set DOPPHWEM to line spacing
             else
               dopphwem=4.3014e-7*sngl(freq)*sqrt(t(i)/molewt(ispeci))
             endif

c  pbhw = abhw*(1-vmr) + sbhw*vmr
             if(kgas.eq.7) then
c             write(*,*) 'Heree',i,vmr(ispeci,i),vmr_tot(kgas,i)
                pbhw=abhw+(vmr_tot(kgas,i)-0.21)*(sbhw-abhw)/(1-0.21)
             elseif(kgas.eq.41) then
                pbhw=abhw+(vmr_tot(kgas,i)-0.79)*(sbhw-abhw)/(1-0.79)
             else
                pbhw=abhw+vmr_tot(kgas,i)*(sbhw-abhw)
             endif

c             if(kgas.eq.7) then
c                pbhw=abhw+(vmr(ispeci,i)-0.21)*(sbhw-abhw)/(1-0.21)
c             elseif(kgas.eq.41) then
c                pbhw=abhw+(vmr(ispeci,i)-0.79)*(sbhw-abhw)/(1-0.79)
c             else
c                pbhw=abhw+vmr(ispeci,i)*(sbhw-abhw)
c             endif

c  Add H2O-broadening, assumed to be 1.35 times the air-broadening.
c  This 1.35 factor has been determined empirically by minimizing
c  the sensitivity of TCCOB xCO2 to retrieved H2O.
c  According to Hartmann, 1.80 would be more realistic.
             if(kgas.ne.1 .and. kgas.ne.49) pbhw=pbhw+
     &         0.35*abhw*vmr(specindex(1)+1,i)

c  Add CO2-broadening, assumed to be 1.15 times the air-broadening
c  This is only important for Venus and Mars
             if(kgas.ne.2) pbhw=pbhw+0.15*abhw*vmr(specindex(2)+1,i)

c Now compute the absorption cross-section of this particular line.
             y=pbhw*p(i)*trat(i)**tdpbhw/dopphwem  ! Pressure broadending (in DW)
             sxcgs=stren*d(i)*vmr(ispeci,i)*vpf(ispeci,i)
             if(scia) sxcgs=sxcgs*d(i)*vmr(ispeci,i)
             if(fcia) sxcgs=sxcgs*d(i)*(1-vmr(ispeci,i))
             eptf=eprime*tfac(i)
             if(eptf.gt.85.0) eptf=85.0
             sxcgsorpidw=sxcgs*exp(eptf)/dopphwem/srpi
c             write(*,*)'sg=',i,sxcgs,eptf,exp(eptf),dopphwem
             godw=grid/dopphwem   ! primative grid point spacing (Doppler widths)
             flinwid=dble((3+amin1(y/frac,
     &       sqrt(abs(y*sxcgsorpidw/srpi/tnulst))))/godw)
c             if(kgas.eq.2) flinwid=dmin1(flinwid,21950.d0) ! Sub-Lorentzian CO2
             kv1=max0(1,1+int(vcent-flinwid)) ! Start index (primative grid points)
             kv2=min0(int(vcent+flinwid),ncp) ! End index (primative grid points)
             nv=kv2-kv1+1    ! Number of primative grid points encompassed by line
             if(nv.gt.0) then
                x1=godw*(dble(kv1)-vcent)  ! starting point (doppler widths)
                if(kgas.eq.2) then
                     sublor=10.5/dopphwem
                else
                     sublor=100.0/dopphwem
                endif
                call gct2_humlik(nv,x1,godw,y,sublor,vvoigt(kv1))
                do jv=kv1,kv2  ! Innermost loop over line profile
                   vac(jv,i,mm)=vac(jv,i,mm)-sxcgsorpidw*vvoigt(jv)
                end do
             endif   ! nv.gt.0
             endif   ! p(i).lt.1.1
          end do     ! i=1,nlev  Loop over levels
        end do       !  do jline=kline1+1,kline2
        close(lun_ll)
      end do         ! end of loop on line files
      write(lunw_rpt,*)nltnu,' lines (all linelists) with S > TNULST'
      close(lunw_rpt)

c  Replace internally-calculated absorption coefficients by externally calculated ones.
      lr=lnbc(runlab)
      lmfile=''
      do jtarg=1,mmax
      if(targmol(specindex(2)+1).eq.jtarg) then
         lmfile=runlab(:10)//'_co2'
         exit
      endif
      if(targmol(specindex(6)+1).eq.jtarg) then
         lmfile=runlab(:10)//'_ch4'
         exit 
      endif
      if(targmol(specindex(7)+1).eq.jtarg)then
         lmfile=runlab(:10)//'_o2'
         exit
      endif
      end do
c      write(*,'(a)') lmfile
c      write(*,*) lnbc(lmfile)
      write(lmfile(lnbc(lmfile)+1:),'(a1,i5.5,a1,i5.5,a7)')
     & '_',nint(fzero+grid),'_',nint(fzero+ncp*grid),'_lm.bin'
c      write(*,*)'lmfile=',lmfile
      inquire(file=lmfile,exist=lmfile_exist)
      if(lmfile_exist) then
c         write(*,*) ' Calling FSIB:'//lmfile
         fsib=file_size_in_bytes(lun_rbin,lmfile)
         write(*,*)' lmfile found: ',lmfile(:lnbc(lmfile))
c         if(fsib.ne.4*ncp*nlev) then
         if(fsib.ne.12*ncp*nlev) then
             write(*,*) lmfile
             write(*,*) 'ncp, nlev, fsib=',ncp,nlev,fsib
             stop 'fsib .ne. 4*ncp*nlev '
         endif
         open(lun_rbin,file=lmfile,access='direct',status='old',
     &   recl=12,form='unformatted')
         irec=0
         do ilev=1,nlev
            do j=1,ncp
               irec=irec+1
               read(lun_rbin,rec=irec)vv,yy,zz
c               vac(j,ilev,jtarg)=-vv-yy-zz                ! CH4: Pine Voigt + LM 
c               vac(j,ilev,jtarg)=vac(j,ilev,jtarg)-yy-zz  ! O2: Hartnam LM & CIA 
               vac(j,ilev,jtarg)=vac(j,ilev,jtarg)-yy     ! O2: Hartnam LM 
c               vac(j,ilev,jtarg)=vac(j,ilev,jtarg)-yy     ! first order LM only
            end do
         end do
         close(lun_rbin)
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
c      open(lun_vac,file=vacfile,status='unknown')
c      write(lun_vac,*)2,nlev
c      write(lun_vac,'(a6,80i12)')' Freq ',(ilev-1,ilev=2,nlev)
c      write(*,*) jtg,ispeci,nlev,ncp
c      do jcp=1,ncp   !  Loop over frequencies
c         write(lun_vac,'(f14.6,8e12.5)')fzero+jcp*grid,
c     &   (-vac(jcp,ilev,jtg)/vmr(ispeci,ilev)/d(ilev),ilev=2,nlev)
c      end do
c      close(lun_vac)
c      end do   !  Loop over target gases

      return
      end

      function vibpf(temp,vibfrq,dgen,nvmode)
c Calculates the vibrational partition function:
c the fraction of molecules that are in the vibrational ground state.
      integer*4 imode,nvmode,dgen(nvmode)
      real*4 temp,vibfrq(nvmode),vibpf
      vibpf=1.0
      do imode=1,nvmode
        vibpf=vibpf*(1.-exp(-(1.43881*vibfrq(imode)/temp)))**dgen(imode)
      end do
      return
      end

