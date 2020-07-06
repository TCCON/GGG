      subroutine spectrum_loop(winfo,debug,
     & lunw_col,lunw_cbf,lcl,colabel,colfile_format,lspmax,
     & runlog,akpath,rayfile,mavfile,targmol,tll_file,
     & parfile,
     & apx,apu,dplist,iptg,ipcl,ipfs,ipsg,ipzo,ipcf,
     & ntg,ncbf,nfp,
     & speci,nspeci_iso,solarll,pars,sptfile)
c
c   Inputs:
c     ntg           I*4  number of target gases
c     speci(ntg)    I*4  The species # for the target gases/isotopologs
c     pars(ntg)     C**  The names of the target gases/isotopologs
c     nspeci_iso    I*4  Number of species in isotopologs.dat
c     nmp           I*4  Number of measured points (in spectrum)
c     ncbf,         I*4  number of continuum terms (basis functions)
c     nfp,          I*4  number of fited parameters = ntg+ncbf+n
c     obsrvd(nmp)   R*4  Measured spectrum (y)
c     apx(nfp)      R*4  A priori state vector (xa)
c     apu(nfp)      R*4  A priori state vector uncertainties (SQRT(diag(Sa)))
c
c   Outputs:
c     lspmax        I*4  Length of longest spectrum named
c     calcul(nmp)   R*4  Calculated spectrum (f(x))
c     cx(nfp)       R*4  State vector (x)
c     ex(nfp)       R*4  State vector uncertainties

c  mode=0    Checks that all input files are readable before investing
c            a lot of time on doing the VAC calculation.
c  mode=1    Does full calculation
c
c  Synthetic spectra are computed when mit=0 and the spectra
c  in the runlog don't exist.
c
c  Pseudo-Code Summary of spectrum_loop.f
c      open(lunr_mav,mavfile)
c      reaf(lunr_mav,specname_mav_
c      oaflag=1
c      mavflag=0
c      ncall=0
c      nspectra=0
c      do ispec=1,mspec
c         read_runlog(data_record(specname_rl)
c         read rayfile(specname_ray)
c         if(specname_rl .ne. specname_ray) stop 'mismatch'
c         kspflag=0
c         call gindfile(dplist,specname_rl,specpath)
c         if(lnbc(specpath).eq.0) kspflag=2
c         hwils=nscycle*resn 
c         if(kspflag.eq.2) hwils=0.0d0
c         if(window exceeds disk file limits) kspflag=1
c         if(specname_rl.eq.specname_mav) then
c            read_mavfile_body(specname_mav)
c            oaflag=1
c            mavflag=1
c         endif
c         if(mavflag.le.0) stop 'Mismatch: specname_rl & specname_mav'
c
c         if(kspflag.ne.1 .and. oaflag.ge.1) then
c            if(ncall.ge.1) then
c               call abscoj()
c            endif          ! (ncall.ge.1) then
c            mav_count=mav_count+1
c            oaflag=0   ! Absorption coefficient in memory are up to date
c         endif          ! (kspflag.eq.0 .and. oaflag.eq.1) 
c
c         if((mit.gt.0.and.kspflag.gt.0) .or.
c     &      (mit.eq.0.and.kspflag.eq.1)) cycle
c
c         if(ncall.ge.1) then
c            if(kspflag.eq.0) then
c               call jetspe(specpath,startm,gint,nmp)
C               if(nmp.lt.ncbf) cycle
c               kspflag=1
c            endif
c
c            if(kspflag.eq.2) then
c               startm=nus
c               gint=graw/interp
c               nmp=1+int((nue-nus)/gint)
c            endif
c
c            call compute_ils
c            call do_retrieval
c            call write_spt
c            call write_col
c         endif ! ncall.ge.1
c         nspectra=nspectra+1
c      end do   !  ispec=1,mspectra 
c      ncall=ncall+1

      implicit none
      real*4 zero, unity, spi

      parameter(zero=0.0)
      parameter(unity=1.0)
      parameter(spi=3.14159265)

      include "ggg_int_params.f"
      include "const_params.f"
      include "int_params.f"

      logical  debug,calc_lm,calc_nv

      integer*4 lf,fbc,reclen_solarll,lso,flagso,idum,i,ii,lsos,jt,
     & ncbf,            ! number of continuum terms (basis functions)
     & imode,           ! 0 = processing mode (do not stop)
     & lsn,lsnd,
     & nfov,kfov,
     & ncall,           ! counts the number of times that subroutine is called
     & oaflag,          ! Indicates whether absorbtion coefficients currently
c                         in memory are obsolete (=1) or not (=0).
     & lcolon,
     & totit,nn,ncol,nspectra,nspectra0,mspectra,
     & ifcsp,ifmsp,iptg,ipcl,ipfs,ipsg,ipzo,ipcf,
     & mcp,ncp,
     & mva,nva,         ! The number of precomputed anscoption coefficients
     & istat,ifm,
     & nscycle,
     & mspt,ispec,
     & freq_flag,       ! =1  presents spectral fits in solar rest frame. 
c                       ! =0  presents spectral fits in atmosphere rest frame.
     & lunw_ak,
     & lun_sts,         ! Solar Transmittance Spectrum
     & lnbc,
     & nspeci_iso,jspeci,
     & nhw,ldec,nmpfp,
     & nmp,imp,
     & mii,            ! Dimension of interpolated (by factor LDEC) ILS
     & ntg,jtg,
     & nfp,ifp,
     & jsp,jva,
     & lcl,lsp,lspmax,
     & interp,
     & defapo,apo_c,apo_m,rc,j,
     & mspxv,
     & mslpd,
     & iseed,
     & mavflag,
     & nlev_ray,ilev,
     & nlev_mav,
     & mit,nit,
     & iyr,iset,
     & kspflag,ifirst,ilast,bytepw,possp,
     & kcp1,kcp2,
     & nsh,nhwmax,
     & lunw_cbf,lunw_col,lunw_spt,lunr_mav,lunr_ray,lun_rlg,mav_count

      parameter (mva=300000000,mcp=3100000,
     & nscycle=25,
     & mii=153847,
     & mslpd=10*mmp,
     & mspxv=14*mcp)

      parameter (lun_sts=24,lun_rlg=25,lunw_ak=26,
     & lunr_ray=27,lunr_mav=28,lunw_spt=29)

      integer*4
     & speci(ntg),
     & targmol(nspeci_iso)

      real*4
     & cfamp,cfperiod,cfphase,
     & vac(mva),
     & ssnmp(mmp),
     & slit(mii),cx(nfp),ex(nfp),
     & obsrvd(mmp),calcul(mmp),cont(mmp),
     & overcol(ntg),oloscol(ntg),
     & zmin,zminwas,
     & fsn,ffr,
     & gasdev,  ! random number generator
     & apx(nfp),apu(nfp),
     & aprx(nfp),apru(nfp),
     & ynoise,
     & rmsocl, 
     & peff,
     & cont_level,cont_tilt,cont_curv,xzo,
     & sza_ray,sza_raywas,bend,
     & tot,
     & vmr,
     & tcp,tcv,
     & wmf(mspeci*mlev),   ! Wet Mole Fraction (aka VMR)
     & vpf(mspeci*mlev),
     & spver(mlev),splos(mlev),
     & cp(mlev),
     & z(mlev),t(mlev),p(mlev),d(mlev),
     & corrld,           ! Factor=[apodixed resolution]/[spectral point spacing
c     & solar_gas_shift,
c     & fovcf,
     & tbar, rnoise,detnoise,sphnoise,
     & slpd(mslpd),
     & pd((mmp+mfp)*mfp), spts(mcp),
     & spxv(mspxv),
     & solzen,roc,fbar,fhw,sfnois,
     & xfs,xsg,rdum

      real*8
     & wlimit,
     & temp,sigma,ptot,ph2o,pn2,ctn2_overtone,ctn2_fundamental,
     & dopp,xCH4,Abs_HT(mcp),
     & riair,
     & hwils,
     & tottc(mtg),tottc2(mtg),toterr(mtg),rmspc,trmspc,wt,twt,
     & frac,             ! fractional size of FOVO compared with solar diameter
     & effres,           ! real spectral resolution (cm-1)
     & resn,             ! 0.5d0/opd = half width of SINC function in cm-1
     & resmax,           ! maximum allowed value of resn
     & opdmin,           ! minimum allowed value of opd
     & opdmax,           ! maximum allowed value of opd
     & rect,             ! frqcen*(fovi**2+amal**2)/8 = width of rectangle(cm-1)
     & resnog,           ! RESN / GRID = 0.5/(OPD*grid)
     & rectog,           ! RECT / GRID
     & rdec,             ! Ratio: GINT/GRID = spectral/primitive point spacings
     & frqcen,           ! centRal frequency (cm-1) of spectral window
     & width,            ! width (cm-1) of spectral window
     & startm,           ! frequency of first returned measured point
     & nus,              ! microwindow starting frequency (cm-1)
     & nue,              ! microwindow ending frequency (cm-1)
     & graw,             ! spacing of raw spectrum (cm-1) from GETINFO
     & gint,             ! spacing of OBSRVD (cm-1) after interpolation
     & grid              ! spacing of primitive spectrum 

      real*8 
     & oblat,            ! observation latitude (deg).
     & oblon,            ! observation longitude (deg).
     & obalt,            ! observation altitude (km)
     & zpdtim,           ! Time of ZPD (UT hours)
     & asza,             ! astronomical solar zenith angle (unrefracted)
     & zenoff,           ! zenith pointing offset
     & apz,              ! asza+zenoff
     & azim,             ! azimuth angle
     & osds,             ! Observer-Sun Doppler Stretch (ppm)
     & fovi,             ! Internal angular diameter of FOV (radians)
     & fovo,             ! External angular diameter of FOV (radians)
     & amal,             ! angular misalignment of interferometer (radians)
     & zoff,             ! zero-level offset
     & snr,              ! Signal-to-Noise Ratio
     & tins,             ! Temperature INSide the instrument
     & pins,             ! Pressure INSide the instrument
     & hins,             ! Humidity INSide the instrument
     & tout,             ! Temperature OUTside the instrument
     & pout,             ! Pressure OUTside the instrument
     & hout,             ! Humidity OUTside the instrument
     & sia,              ! Solar Intensity Average (arbitrary units)
     & fvsi,             ! Fractional Variation in Solar Intensity
     & wspd,             ! Wind Speed
     & wdir,             ! Wind Direction
     & aipl,             ! Airmass-Independent Path Length (km)
     & lasf,             ! laser frequency (e.g. 15798.03 cm-1)
     & wavtkr,           ! suntracker operating frequency (e.g. 9900 cm-1)
     & opd               ! Optical path difference (cm) of interferogram

      character winfo*(*),specname_rl*(nchar),pars(ntg)*(*),
     & data_fmt_read_rl*256,
     & col_labels_rl*320,
c     & ss(nfp)*4,
     & sptfile*(*),akpath*(*),specpath*(mfilepath),
     & sptpath*(mfilepath),header*256,rayfile_format*256,
     & solarll*(*),specname_mav*(nchar),
     & path_n2_cia_fundamental*128, path_h2o_n2_cia_fundamental*128,
     &  path_n2_cia_overtone*128,
     & string*80,colfile_format*(*),cbffile_format*40,
     & oformat*12,colabel*(*),
     & mavstring*(nchar+14),tll_file*(*),parfile*(*),dplist*(*),
     & col1*1,apf*2,rayfile*(*),specname_ray*(nchar),
     & runlog*(*),mavfile*(*)

      parameter (mspectra=2048*1024,resmax=0.625d0,opdmin=0.5d0/resmax)

      save ispec,ncall,nspectra0
      data ispec/0/
      data ncall/0/

      rdum=big      ! Prevent compiler warning (unused variable)
      idum=mauxcol  ! Prevent compiler warning (unused variable)
      idum=mcolvav  ! Prevent compiler warning (unused variable)
      idum=mgas     ! Prevent compiler warning (unused variable)
      idum=mrow_qc  ! Prevent compiler warning (unused variable)
      idum=mtg      ! Prevent compiler warning (unused variable)
      idum=mvmode   ! Prevent compiler warning (unused variable)
      rdum=tiny     ! Prevent compiler warning (unused variable)

      rayfile_format='(a,200f11.0)'  ! to support old .ray files.

c      path_n2_cia_fundamental='/home/toon/ddd/linelist/CT-N2.N2'
c      path_h2o_n2_cia_fundamental='/home/toon/ddd/linelist/CT-N2.H2O'
c      path_n2_cia_overtone=
c     & '/home/toon/ddd/linelist/N2_in_air_CIA_0-2_band.dat'
      calc_lm=.false.
      oaflag=1
      mavflag=0
      imode=0
c  Read max # of SPT files (if a value is provided on the SPT line of the .ggg file)
      mspt=4000  ! default value
      lf=fbc(sptfile)
      if(lnbc(sptfile).gt.lf) read(sptfile(lf:),*) mspt

c      lc=index(winfo,':')
c      call substr(winfo(lc+1:),pars,mtg,ntg)
c      if(ntg.gt.mtg) then
c          write(*,*)' spectrum_loop: Error: NTG > MTG ',ntg,mtg
c          stop 'Increase parameter MTG inspectrum_loop.f '
c      endif

c      if( index(winfo,'debug') .gt. 0 ) then
c         debug=.true.
c      else
c         debug=.false.
c      endif

      if (debug) write(*,*) 'Entered spectrum_loop: ncall=',ncall

      if(ncall.ge.1) then
         write(oformat,'(a5,i2.2,a5)')'(1x,a',lspmax+1,',a,a)'
c         write(*,*)'oformat = ',oformat
         write(lunw_col,oformat)' Spectrum            ',
     &   colabel(:lnbc(colabel))
      endif 

      open(lunr_mav,file=mavfile,status='old')
      read(lunr_mav,*)
      read(lunr_mav,'(14x,a)')string  ! Next Spectrum
      lcolon=index(string,':')   ! Next Spectrum:
      read(string(lcolon+1:),'(a)') specname_mav   ! First Next Spectrum
c      write(*,*)'Next Spectrum = ',specname_mav

      open(lunr_ray,file=rayfile,status='old')
      read(lunr_ray,*)nn,ncol
      do j=2,nn
         read(lunr_ray,'(a)') header
         if(header(:7).eq.'format=') rayfile_format=header(8:)
      end do
c      call skiprec(lunr_ray,nn-1)
      nlev_ray=ncol-7  ! First 7 columns are (spec,Zobs,Pobs,SZA,Bend,FOV,Zmin)
c
      if(nlev_ray.gt.mlev) then
         write(*,*) 'nlev_ray,mlev=',nlev_ray,mlev
         stop 'Increase parameter MLEV_ray'
      endif
 
      read(winfo,*) frqcen,width,mit,defapo,interp,freq_flag
c  INTERP is factor by which the observed spectrum will be over-sampled
c  LDEC is factor by which ILS is oversampled for improved accuracy.
      grid=0.666666d-06*frqcen
      nus=frqcen-width/2
      nue=frqcen+width/2
12    kcp1=int(nus/grid)+1
      kcp2=int(nue/grid)
      nsh=int(2+2*resmax/grid)
      nhwmax=nint(nscycle*resmax/grid)
      ncp=kcp2-kcp1+1+2*nhwmax+2*nsh
      ifcsp=kcp1-nhwmax-nsh  ! Index of First Calculated Spectral Point
c      write(*,*) 'nus, nue, grid, ncp, nhwmax, nsh= ',
c     & nus, nue, grid, ncp, nhwmax, nsh
      nva=ncp*nlev_ray*(ntg+1)
      if(nva+ncp.gt.mva) then
         write(6,'(2(a,i5),a)')'Increase MVA from',mva/1000000,
     &   'M  to',(nva+ncp)/1000000,'M'
         write(6,*)'to avoid loss of accuracy'
         grid=1.0022d0*grid*(nva+ncp)/mva
         go to 12
      endif

      if(ncp.gt.mcp) then
         write(*,*) 'MCP, NCP=',mcp,ncp
         stop 'Increase parameter MCP'
      endif

      if(nfp*ncp.gt.mspxv) then
         write(6,*)' Increase MSPXV=',mspxv,' to ',ncp*nfp
         stop 'spectrum_loop: Increase parameter MSPXV'
      endif
C
      do jt=1,ntg
         tottc(jt)=0.0d0
         tottc2(jt)=0.0d0
         toterr(jt)=0.1d-36
      end do
      twt=0.0d0
      trmspc=0.0d0
      totit=0
      mav_count=0
c      write(*,*) 'runlog=',runlog
      open(lun_rlg,file=runlog,status='unknown')
      call read_runlog_header(lun_rlg,data_fmt_read_rl,col_labels_rl)
c      read(lun_rlg,*,err=888) nlhead,ncol
c      do i=2,nlhead
c         read(lun_rlg,*)
c      end do
c888   continue  !  Continue to support old format runlogs
      nspectra=0
      if(debug) write(*,*)' Main loop...',nspectra,mspectra
      do ispec=1,mspectra         !  Main fitting loop over spectra
         call read_runlog_data_record(lun_rlg,data_fmt_read_rl,
     &   col1,specname_rl,iyr,iset,zpdtim,
     &   oblat,oblon,obalt,asza,zenoff,azim,osds,
     &   opd,fovi,fovo,amal,ifirst,ilast,graw,possp,bytepw,zoff,snr,apf,
     &   tins,pins,hins,tout,pout,hout,
     &   sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)
c         write(*,*) 'ispec,specname_rl =',ispec,specname_rl
         if(debug) write(*,*) 'specname_rl, istat=',specname_rl,istat
         if(istat.eq.3) write(*,*)
     &   'read_runlog_data_record: format err: istat=',istat,specname_rl
         if(istat.ne.0) exit
         if(col1.eq.':') cycle

         opdmax=dabs(0.5/graw)
         lsp=lnbc(specname_rl)
         specname_rl=specname_rl(:lsp)
c
         read(lunr_ray,rayfile_format) specname_ray,rdum,rdum,sza_ray,
     &   bend,rdum,zmin,(splos(j),j=1,nlev_ray)
c         write(*,*) 'specname_ray = ',specname_ray
c         write(*,*) 'zmin=',zmin
c         write(*,*)'nlev_ray,splos=',nlev_ray,specname_ray,
c     &  (splos(j),j=1,nlev_ray)
c         write(37,*)zmin,zminwas,sza_ray,sza_raywas,
c     &   (zmin-zminwas)/(sza_ray-sza_raywas)
         zminwas=zmin
         sza_raywas=sza_ray
         specname_ray=specname_ray(:lsp)
         if(specname_ray.ne.specname_rl) then
            write(6,'(a,2x,a)') ' specname_ray        specname_rl'
            write(6,'(a,2x,a)') specname_ray,specname_rl
            stop 'spectrum mismatch 2'
         endif

c  If MIT > 0, check that requested spectrum is on disk,
c  and that it covers the specified spectral interval.
         kspflag=0
         call gindfile(dplist,specname_rl,specpath)
c         write(*,*) 'After gindfile: ',specname_rl,specpath
         if(ispec.eq.1 .and. ncall.eq.0) write(*,*)'specpath = '//
     &   specpath(:lnbc(specpath))
         if(lnbc(specpath).eq.0) then
            kspflag=2
c            write(*,*) 'Spectrum not found: '//specname_rl
         endif
         if(debug) write(*,*)' kspflag = ',
     &     kspflag,'  '//specname_rl,specpath
c         write(*,*)'H kspflag = ',kspflag,'  '//specpath(:32)
c
c  Apply air-to-vacuum correction
c         write(*,*)'tins,pins,hins',tins,pins,hins
         if(kspflag.lt.2) graw=graw*riair(lasf,tins,pins,hins)/
     &   riair(frqcen,tins,pins,hins)

c  Apply FOV corrections
c         graw=graw*(1.D0+(amal**2+fovi**2)/16)  ! FOV correction

         resn=0.5d0/opd
         if(opd.gt.opdmax) then
            resn=0.5d0/opdmax
            write(*,*)'Warning: OPD > OPDmax = 0.5/graw: ',
     &      specname_rl(:lsp),opd,opdmax
            if (imode.gt.0) stop 'OPD > '
         endif
         if(opd.lt.opdmin) then
            resn=0.5d0/opdmin
            write(*,*)'Warning: OPD < OPDmin = 0.5/RESmax: ',
     &      specname_rl(:lsp),opd,opdmin
            if (imode.gt.0) stop 'OPD < OPDmin'
         endif

c         if(1.0001*resn.lt.graw) then
c            write(*,*) 'resn,graw=',resn,graw,specname_rl
c            stop ' resn < graw'
c         endif

c  Measured spectrum must exceed fitting interval to
c  allow convolution with ILS
         hwils=nscycle*resn ! half-width of the slit function in cm-1 (always +ve)
         fbar=0.5*sngl(graw)*(ifirst+ilast)
         fhw=0.5*abs(sngl(graw))*(ilast-ifirst)
         if(kspflag.eq.2) hwils=0.0d0               ! Making synthetic spectrum
         if(debug)write(*,*)kspflag,hwils,nus,nue,fbar,fhw
         if(debug)write(*,'(a,3f13.6)')'nus, nue, hwils =',nus,nue,hwils
         if(debug)write(*,'(a,2f13.6)')'w_ifirst,w_ilast=',graw*ifirst,
     &   graw*ilast
c         if( nus-hwils .lt. graw*ifirst ) kspflag=1   ! Lower window limit < disk file
c         if( nue+hwils .gt. graw*ilast  ) kspflag=1   ! Upper window limit > disk file
         if( nus-hwils .lt. fbar-fhw ) kspflag=1   ! Lower window limit < disk file
         if( nue+hwils .gt. fbar+fhw ) kspflag=1   ! Upper window limit > disk file
c========================================================================
c  Read model & wmf information (SUNRUN.MAV)
         if(kspflag.ne.1 .and. lsp.gt.lspmax) lspmax=lsp  ! Longest fitted spectrum name
c          write(*,*)'specname_rl, specname_mav = ',
c     &    specname_rl(:24),specname_mav(:24),kspflag,oaflag
         if(specname_rl.eq.specname_mav) then
c           write(*,*)'Calling mavfile body',kspflag,oaflag,mavflag
            call read_mavfile_body(lunr_mav,nlev_ray,nspeci_iso,
     &      nlev_mav,z,t,p,d,wmf,parfile)
c           write(*,*)'Exited mavfile body',kspflag,oaflag,mavflag
c            if (nlev_mav.ne.nlev_ray) write(*,*)'Warning: spectrum_loop:'//
c     &     '  nlev_mav .ne. nlev_ray',nlev_mav,nlev_ray
            mavflag=1
            read(lunr_mav,'(a)',end=66) mavstring
c            write(*,*)'mavstring= ',mavstring,specname_rl
            if(mavstring(1:14).eq.'Next Spectrum:') then
               read(mavstring(15:),'(a)') specname_mav
c               write(*,*)'Next Spectrum = ',specname_mav
            else
               write(*,*) mavstring
               stop 'Failed to find Next Spectrum: key-word'
            endif
66          continue
            if(index(winfo,' sa_temp ').gt.0)
     &      call vadd(t,1,4.,0,t,1,nlev_ray) 
            if(index(winfo,' sa_pres ').gt.0)
     &      call vmul(p,1,.95,0,p,1,nlev_ray)
            oaflag=1  ! Absorption coefficients currently in memory are obsolete
         endif  ! specname_rl.eq.specname_mav

         if(mavflag.le.0) then
            write(*,*)'Mismatch between spectra in mavfile & runlog'
            write(*,*) specname_mav,specname_rl
            stop 'mavfile data have not been read in -- cannot proceed.'
         endif

         if(debug)write(*,*)'ncall,kspflag,oaflag=',ncall,kspflag,oaflag
c         if(kspflag.eq.0 .and. oaflag.ge.1) then  !  current spectrum encompasses window
         if(kspflag.ne.1 .and. oaflag.ge.1) then  !  current spectrum encompasses window
c  Pre-compute absorption coefficient
            if(ncall.ge.1) then  ! Do time-consuming stuff, skipped on previous call
c               write(*,*)' Calling abscoj:',nlev_mav,t(1),p(1)
               if(debug) write(*,*)' Calling abscoj:',nlev_mav,t(1),p(1)
               call vmov(zero,0,vac,1,nva)
               if(index(winfo,' lm ').gt.0) calc_lm=.true.
               if(index(winfo,' nv ').gt.0) calc_nv=.true.
               call abscoj(specname_rl,nlev_mav,t,p,d,nspeci_iso,
     &         targmol,wmf,vpf,tll_file,parfile,ifcsp,grid,ncp,vac,
     &         vac(nva+1),calc_lm,calc_nv)
               write(6,'(a4,'//oformat(2:))' n-i','  Spectrum        ',
     &         colabel(:lcl+40)
               if(debug)  write(*,*)' Called abscoj...'

c  Code to implement Ha Tran's CH4 lineshape
               if(2.eq.1) then
                  write(*,*)' Tran CH4 R6 pCqSDHC: ',ifcsp*grid,
     &            (ifcsp+ncp)*grid
                  if(ifcsp*grid.lt.6078.0 .and.
     &            (ifcsp+ncp)*grid.gt.6076.2) then
c                     ii=1              ! Non-target gas
                     ii=1+ncp*nlev_mav ! 1'st target gas
                     do ilev=1,nlev_mav
                        temp=t(ilev)
                        ptot=p(ilev)
                        xCH4 = wmf(speci(1)+(ilev-1)*nspeci_iso)
                        write(*,*) ilev,temp,ptot,xCH4
                        Call CompAbs(ifcsp*grid,(ifcsp+ncp)*grid,grid,
     &                  temp,ptot,xCH4,Abs_HT)
                        do i=1,ncp
                           write(55,*)ilev,i,(ifcsp+i)*grid,Abs_HT(i)
                           vac(ii)=vac(ii)-sngl(Abs_HT(i))
                           ii=ii+1
                        end do
                     end do ! ilev=1,nlev_mav
                  endif   !  if(ifcsp*grid.lt.5000.0 .and. 
               endif 

c  Code to implement Hartmann's N2 CIA
               if(1.eq.2) then
                  write(*,*) 'Using Hartmanns N2 CIA'
                  if(ifcsp*grid.lt.5000. .and.
     &            (ifcsp+ncp)*grid.gt.4300.) then
                     call lecn2_overtone(path_n2_cia_overtone)
c                     ii=1              ! Non-target gas
                     ii=1+ncp*nlev_mav ! 1'st target gas
                     do ilev=1,nlev_mav
                        temp=t(ilev)
                        ptot=p(ilev)
                        pn2=0.7808*ptot
                        do i=1,ncp
                           sigma=(ifcsp+i-1)*grid
                           vac(ii)=vac(ii)-sngl(ctn2_overtone(sigma,pn2,
     &                     ptot,temp))
                           ii=ii+1
                        end do
                     end do ! ilev=1,nlev_mav
                  endif   !  if(ifcsp*grid.lt.5000.0 .and. 

                  if(ifcsp*grid.lt.2830. .and.
     &            (ifcsp+ncp)*grid.gt.1930.) then
                     call lecn2_fundamental(path_n2_cia_fundamental,
     &               path_h2o_n2_cia_fundamental)
                     ii=1              ! Non-target gas
                     ii=1+ncp*nlev_mav ! 1'st target gas
c                     write(*,*)' WMF ilev=3:',wmf(1+(3-1)*nspeci_iso),p(3),
c     &               0.7808*p(3)*(1-wmf(1+(3-1)*nspeci_iso))
                     do ilev=1,nlev_mav  ! Loop over atmospheric levels
                        temp=t(ilev)
                        ptot=p(ilev)  
                        ph2o=ptot*wmf(1+(ilev-1)*nspeci_iso)
c                        ph2o=0.0  !  Uncomment to zero N2-H2O absorption
                        pn2=0.7808*(ptot-ph2o)
                        do i=1,ncp
                           sigma=(ifcsp+i-1)*grid
                           vac(ii)=vac(ii)-sngl(ctn2_fundamental(sigma,
     &                     pn2,ph2o,ptot,temp))
                           ii=ii+1
                        end do
                     end do ! ilev=1,nlev_mav
                  endif   !  if(ifcsp*grid.lt.2830.0 .and. 

               endif   !  if(1.eq.2) then

            endif          ! (ncall.ge.1) then
            mav_count=mav_count+1
            oaflag=0   ! Absorption coefficient in memory are up to date
         endif          ! (kspflag.eq.0 .and. oaflag.eq.1) 


         if((mit.gt.0.and.kspflag.gt.0) .or.
     &      (mit.eq.0.and.kspflag.eq.1)) then
            if(debug) then
               write(*,*)' Requested interval: ',nus,' to ',nue
               write(*,*)' Required interval:',nus-hwils,' to',nue+hwils
               write(*,'(a,a,f9.3,a4,f9.3,a5)') specname_rl,
     &       ' only encompasses',ifirst*graw,' to ',ilast*graw,' cm-1'
            endif
            cycle ! skip missing/partial spectrum
         endif
c=========================================================
c         resn=0.5d0/opd
c         if(resn.lt.graw) resn=graw
         if(index(winfo,' sa_fovi ').gt.0) fovi=fovi*1.07
         rect=frqcen*(fovi**2+amal**2)/8  ! old code
         if(ncall.ge.1) then
c
c  This code from DG is a crude attempt to simulate
c  a frequency-independent misalignment effect.
c         if(runlog(ir+1:ir+2).eq.'in') then
c            rect=frqcen*(fovi**2)/8 + 2100.D0*(amal**2)/8
c         elseif(runlog(ir+1:ir+2).eq.'hg') then
c            rect=frqcen*(fovi**2)/8 + 800.D0*(amal**2)/8
c         else
c            rect=frqcen*(fovi**2+amal**2)/8
c         endif
c
c  Select apodization function
            if(apf.eq.'BX') then
               apo_m=defapo
               apo_c=defapo
            else  ! if the measured spectra are already apodized
               apo_m=0  ! for perfect representation of synthetic spectra
               if(apf.eq.'N1') then
                  apo_c=1
               elseif(apf.eq.'N2') then
                  apo_c=2
               elseif(apf.eq.'N3') then
                  apo_c=3
               elseif(apf.eq.'TR') then
                  apo_c=4
               else 
                  write(6,'(a21,1x,a2,1x,a)')
     &            specname_rl,apf,' ???  Unknown apodization function'
                  stop
               endif
            endif

c---------------------------------------------------------
c  FIND the spectrum file, return the PATH to the spectrum
            call vmov(zero,0,obsrvd,1,mmp)
            if(kspflag.eq.0) then
               if(debug) write(*,*)'Calling jetspe: ',specpath,kspflag,
     &         resn,graw,ifirst,ilast,possp,bytepw,nus,nue,apo_m,interp
               call jetspe(specpath,resn,graw,ifirst,ilast,possp,bytepw,
     &         nus,nue,apo_m,interp,zero,zero,
     &         obsrvd,mmp,nmp,startm,gint,rc)
               if(nmp.lt.ncbf) then
                  write(*,*)'skipping spectrum '//specname_rl
                  write(*,*)'nmp < ncbf = ',nmp,ncbf
                  write(*,*)'Widen window or reduce NCBF'
                  cycle
               endif
               if(debug) write(*,*)'obsrvd=',obsrvd(1),obsrvd(2)
c              write(*,*)'obsrvd=',obsrvd(1),obsrvd(2),obsrvd(3)
               if(rc.ne.0) then
                  write(6,*)' Error in JETSPE. Spectrum ',specname_rl,rc
                  kspflag=1
               endif
            endif  ! kspflag.eq.0
            if(kspflag.eq.2) then
               startm=nus
               gint=graw/interp
               nmp=1+int((nue-nus)/gint)
            endif
            ifmsp=nint(startm/gint)  ! index of first measured spectral point
            if(nmp.gt.mmp) then
               write(*,*)'nmp,mmp=',nmp,mmp
               stop 'Increase parameter MMP'
            endif

            lsn=index(winfo,' sa_snr ')
            lsnd=index(winfo,' sa_snr/')
c            write(*,*) 'lnn,lsnd= ',lsn,lsnd
            if(lsn.gt.0 .or. lsnd.gt.0) then  ! Add noise
c  Add noise and systematic error to OBSRVD (assumes a continuum level of 1.0)
               call vdot(obsrvd,1,unity,0,tot,nmp)
               tbar=tot/nmp
      
c  sfnoise represent the (inverse of) the spectral response.
c  necessary because we do everything in transmittance space
c  rather than uncalibrated signal space.
               if (frqcen.lt.1880.) then
                  sfnois=sngl(1+((frqcen-1400)/550)**2)  ! HgCd
                  detnoise=0.0008*sfnois  ! MATMOS detector noise
                  sphnoise=0.0012*sfnois  ! MATMOS source photon noise
               else
                  sfnois=sngl(1+((frqcen-2700)/1200)**2) !  InSb
                  detnoise=0.0004*sfnois  ! MATMOS detector noise
                  sphnoise=0.0017*sfnois  ! MATMOS source photon noise
               endif

               if(lsnd.gt.0) then
                  read(winfo(lsnd+8:),*) fsn ! =5.8 for NOMAD and ACS (QQQQQ)
                  detnoise=detnoise/fsn  !   NOMAD detector noise
                  sphnoise=sphnoise/fsn  !   NOMAD source photon noise
               endif
c            write(*,*) fsn,tbar,detnoise,sphnoise
   
c  Add detector noise and source photon noise in quadrature
               rnoise=SQRT(tbar*sphnoise**2+detnoise**2)
c               write(*,*)'tbar,rnoise,SNR =',tbar,rnoise,tbar/rnoise
               do imp=1,nmp
                  obsrvd(imp)=obsrvd(imp)+rnoise*gasdev(iseed) ! random noise
c                  obsrvd(imp)=obsrvd(imp)*(1+0.005*obsrvd(imp))  !  0.5%  detector non-linearity
               end do
            endif  !  Add noise
            iseed=iseed+1113

            if(index(winfo,' sa_cf ').gt.0) then  ! Add channeling
c  Add noise and systematic error to OBSRVD (assumes a continuum level of 1.0)
               call vdot(obsrvd,1,unity,0,tot,nmp)
               tbar=tot/nmp
c  Assume channeling amplitude is 0.05% of the continuum level 
c  Assume two different frequencies (0.41 cm-1 & 1.01 cm-1)
               do imp=1,nmp
                  obsrvd(imp)=obsrvd(imp)+
     &            0.0004*tbar*cos(2*spi*imp*sngl(gint)/0.41)
                  obsrvd(imp)=obsrvd(imp)+
     &            0.0004*tbar*cos(2*spi*imp*sngl(gint)/1.01)
               end do
            endif  !  Add channelling

            nmpfp=nmp+nfp
c  Pre-compute ILS (oversampled by a factor LDEC) and normalize each
c  interleaved component to unity.
            rdec=gint/grid
c         write(*,*) gint,grid,rdec
            nhw=iabs(nint(nscycle*resn/grid))
            ldec=1+int(abs(16*grid/resn))
            if(1+2*nhw*ldec.gt.mii) stop ' Increase MII'
            resnog=ldec*resn/grid
            rectog=ldec*rect/grid
c            write(*,*)'interp,ldec,nhw=',interp,ldec,nhw
            call compute_ils(apo_c,nhw,ldec,resnog,rectog,0.0d0,slit)

c Write out ILS
c            open(56, file='fort.56')
c            write(56,*) 2,2
c            write(56,*)' v s'
c            do j=1,1+2*nhw*ldec
c               write(56,*)grid*(j-float(1+nhw*ldec))/ldec,slit(j)
c            end do
c            close(56)
c===================================================================
c  Error estimation naively assumes spectral points are linearly independent.
c  CORRLD is an estimate of how dependent neighbouring points really are.
c  CORRLD = 1 means that points are truly linearly independent.
c  Variances must be scaled by CORRLD to ensure that when the same spectrum is
c  analyzed under different values of APO or INTERP, the same size errors result
c            corrld=sqrt(rect**2+((1.+0.4*apo_c)*resn)**2+(resn-0.5d0/opd)**2)
c     &    /gint
            effres=sqrt(rect**2+((1.+0.4*apo_c)*resn)**2)
            corrld=sngl(effres/gint)

c========================================================
c  Compute SPXV
c  The first column of SPXV contains the non-target VACS
c  The next NTG columns of SPXV contain the target VACs
c  The NTG+2'nd column contains the Solar Transmittance Spectrum
c  The NTG+3'rd and NTG+4'rd columns are workspace.
c  SPXV does not change for a given observation geometry,
c  even though the absorber amounts are being iterated.
            jva=1
            jsp=1
            do jtg=0,ntg
               call vmov(zero,0,spxv(jsp),1,ncp)
               do ilev=1,nlev_mav
                  call vsma(vac(jva),1,ckm2cm*splos(ilev),spxv(jsp),1,
     &            spxv(jsp),1,ncp)
c                  write(*,*)'spxv: xxx:',ilev,z(ilev),jva,splos(ilev),vac(jva)
                  jva=jva+ncp
               end do
               jsp=jsp+ncp
            end do  ! jtg=0,ntg

c  Compute effective pressure of first target gas
            tcv=0.0001
            tcp=0.0
            if(ntg.gt.0) then
               do ilev=1,nlev_mav
                  vmr = wmf(speci(1)+(ilev-1)*nspeci_iso)
                  tcv=tcv+d(ilev)*vmr
                  tcp=tcp+d(ilev)*vmr*p(ilev)
               end do
            end if
            peff=tcp/tcv

c   Compute Solar Pseudo-Transmittance Spectrum
c  (if a solar linelist is included in the .ggg.file)
            if(index(winfo,' so ').gt.0.or.index(winfo,' sg ').gt.0 .or.
     &     index(winfo,' so/').gt.0 .or. index(winfo,' sg/').gt.0 ) then
               lso=lnbc(solarll)
               if( lso .gt. 1 ) then
                  flagso=1
                  read(solarll(lso-2:lso),*) reclen_solarll
                  if(reclen_solarll.eq.108) then
                     frac=fovo/9.2e-3  ! the sun is 9.2 mrad in diameter on average
                     call solar_pts(lun_sts,solarll,
     &               ifcsp,grid*(1.d0+1.d-6*osds),frac,spts,ncp)
                  else
                     write(*,*)'Solarll= ',solarll
                     write(*,*)'lso= ',lso
                     write(*,*) 'reclen_solarll= ',reclen_solarll
                     stop 'Upgrade solar linelist to .108 format'
                  endif
               else
                  stop 'No named solar linelist in .ggg file'
               endif
cZZ Restored by ZCZ 12/23/2017 based on GGG2014
               lsos=index(winfo,' so/')
               if(lsos.gt.0) then
                  read(winfo(lsos+4:),*) ffr  ! absorption reduction factor
c                  write(*,*)'Reducing solar absorption depths by factor',ffr
                  do i=1,ncp
                     spts(i)=1.0-(1.0-spts(i))/ffr
                  end do
               endif
cZZ END 
            else
               flagso=0
               osds=0.0d0
               call vmov(unity,0,spts,1,ncp)
            endif  ! lso .gt. 0
            if(index(winfo,' sa_zoff ').gt.0) zoff=zoff+0.005

c===================================================================
c  Compute the vertical slant path distances
c            write(*,*)'nlev_mav,nlev_ray,ncell=',nlev_mav,nlev_ray,ncell
            if (nlev_ray.le.9) then
               do j=1,nlev_ray
                  spver(j)=splos(j)
               end do
            else
               call compute_vertical_paths(ncell,zmin,z,d,spver,
     &         nlev_mav)
            endif
c            do j=1,nlev_ray
c               write(*,*) j,z(j),d(j),spver(j),splos(j)
c            end do
c            write(*,*)'splos=',(splos(j),j=1,nlev_ray)
c            write(*,*)'spver=',(spver(j),j=1,nlev_ray)
c            write(*,*)'  d  =',(d(j),j=1,nlev_ray)
c            write(*,*)' wmf =',(wmf(speci(1)+(j-1)*nspeci_iso),j=1,nlev_ray)
c===================================================================
c Compute vertical and LOS slant columns for the target gases
c Do in reverse order so that the concentration profile of the
c first target gas is retained in CP at the end.
            do jtg=ntg,1,-1
               jspeci=speci(jtg)
               call vmul(wmf(jspeci),nspeci_iso,d,1,cp,1,nlev_mav)
               call vdot(cp,1,spver,1,overcol(jtg),nlev_mav)
               call vdot(cp,1,splos,1,oloscol(jtg),nlev_mav)
               overcol(jtg)=overcol(jtg)*ckm2cm
               oloscol(jtg)=oloscol(jtg)*ckm2cm
c            write(*,*)'jtg,jspeci,overcol,oloscol',jtg,jspeci,
c     &       overcol(jtg),oloscol(jtg)
            end do
c            do i=1,nlev_ray
c              write(*,*)'wmf=',i,z(i),(wmf(speci(jtg)+nspeci_iso*(i-1)),jtg=1,ntg)
c            end do
c
c=====================================================================
c  The structure of the state vector is as follows. Elements:
c      1 to NTG            are the Target Gases
c  NTG+1 to NTG+NCBF       are the Continuum Basis Function Coefficients
c  NTG+NCBF+NFS            is the FS (if NFS>0)
c  NTG+NCBF+NFS+NZO        is the ZO (if NZO>0)
c  
c  Read the initial values of the parameters to be fitted and their estimated
c  a priori variances.
c            open(lun_apx,file=ap_file,status='old')
c            read(lun_apx,*)
c            kcbf=ntg+1
c            do j=1,5       !  Read a priori values & uncertainties
c               read(lun_apx,'(f5.0,f9.0,1x,a4)') apx(kcbf),apu(kcbf),ss(kcbf)
c               if(index(winfo,ss(kcbf)).gt.0) kcbf=kcbf+1
c            end do
c            do jj=ntg+1,ntg+5
c               read(lun_apx,'(f5.0,f9.0,1x,a4)') apx(jj),apu(jj),ss(jj) 
c                   if(ss(jj).eq.' cl ') then
c                   ipcl=jj
c               elseif(ss(jj).eq.' ct ') then
c                   ipct=jj
c               elseif(ss(jj).eq.' cc ') then
c                   ipcc=jj
c               elseif(ss(jj).eq.' fs ') then
c                   ipfs=jj
c               elseif(ss(jj).eq.' zo ') then
c                   ipzo=jj
c               else
c                   write(*,*) ntg,jj,jj-ntg,ss(jj)
c                   stop 'Unknown/Unreadable mnemonic found in a priori file'
c               endif
c            end do

c            read(lun_apx,*) apx(ipcl),apu(ipcl),ss(ipcl) ! Continuum Level
c            read(lun_apx,*) apx(ipct),apu(ipct),ss(ipct) ! Continuum Tilt
c            read(lun_apx,*) apx(ipcc),apu(ipcc),ss(ipcc) ! Continuum Curvature
c            read(lun_apx,*) apx(ipfs),apu(ipfs),ss(ipfs) ! Frequency Shift
c            read(lun_apx,*) apx(ipzo),apu(ipzo),ss(ipzo) ! Zero Offset
c            read(lun_apx,*)                              ! Solar Scaling
c            if(ntg.ge.1) read(lun_apx,*) apx(1),apu(1)   ! First Target Gas
c            if(ntg.ge.2) read(lun_apx,*) apx(2),apu(2)   ! Other Target Gases
c            close(lun_apx)
c            call vmov(apx(2),0,apx(3),1,ntg-2) ! absorber amount
c            apx(ipzo)=apx(ipzo)+sngl(zoff)
c            call vmov(apu(2),0,apu(3),1,ntg-2)  ! wmf of non-target gases
c            apu(ipfs)=apu(ipfs)*sngl(1.d0+rdec) ! shift   up from 0.5 4-DEC-98 gct
c==========================================================================
c            write(*,*)' ntg=',ntg
c            write(*,*)'sl: apx(ipct)=',ipct,apx(ipct)
            call vmov(apx,1,aprx,1,nfp) ! Initialize to A PRIORI values each spectrum
            if(ipzo.gt.0) then
               aprx(ipzo)=aprx(ipzo)+sngl(zoff)
               xzo=aprx(ipzo)
            else
               xzo=sngl(zoff)
            endif
            call vmov(aprx,1,cx,1,nfp)  ! Initialize to A PRIORI values each spectrum
            call vmov(apu,1,apru,1,nfp) ! Initialize to A PRIORI uncertainties each spectrum
            if(ipfs.gt.0) apru(ipfs)=apru(ipfs)+sngl(1.0d0+rdec)
c            if(ipsg.gt.0) apru(ipsg)=apru(ipsg)+sngl(1.0d0+rdec)
ccc            call vmov(aprx,1,cx,1,nfp)  ! Initialize to A PRIORI uncertainties each spectrum
            if(debug)
     &     write(6,*)  'It    CL       CT    CC    FS    ZOFF  RMS/CL'//
     &     '  Vfact1  Vfact2  Vfact3  Vfact4  Vfact5  Vfact6'

            if(index(winfo,' prof_ret_only ').gt.0) ifm=4  ! Profile Retrieval only
            if(index(winfo,' prof_ret+scale ').gt.0) ifm=5  ! Profile Retrieval
            kfov=index(winfo,'nfov=')
            if(kfov.gt.0) then
               read(winfo(kfov+5:),'(i1)')nfov
               ifm=3
               solzen=sngl(asza+zenoff)
               roc=6378.
               fbar=(ifirst+ilast)*sngl(graw)/2
c               write(*,*) 'Calling do_retrieval3...'

               call do_retrieval3(obsrvd,nmp,aprx,apru,slit,nhw,
     &         ifcsp,ifmsp,iptg,ipcl,ipfs,ipsg,ipzo,ipcf,xzo,
     &         cont_level,cont_tilt,cont_curv,
     &         cfamp,cfperiod,cfphase,
     &         nfov,z,t,p,solzen,fovo,roc,obalt,wavtkr,fbar,
     &         ldec,rdec,spts,spxv,vac,splos,nlev_mav,ncp,ntg,ncbf,nfp,
     &         snr,corrld,debug,mit,nit,cont,calcul,rmsocl,
     &         cx,ex,slpd,pd,ssnmp)
            else  !  FOV center only.
               ifm=1
               nfov=1
               if(debug)write(*,*)'Call do_retrievl cx=',(cx(j),j=1,nfp)
               if(debug)write(*,*)'Call do_retrieval spxv(1,0)=',spxv(1)
               call do_retrieval(obsrvd,nmp,aprx,apru,slit,nhw,
     &         ifcsp,ifmsp,iptg,ipcl,ipfs,ipsg,ipzo,ipcf,xzo,
     &         cont_level,cont_tilt,cont_curv,
     &         cfamp,cfperiod,cfphase,
     &         ldec,rdec,spts,spxv,vac,splos,nlev_mav,ncp,ntg,ncbf,nfp,
     &         snr,corrld,debug,mit,nit,cont,calcul,rmsocl,cx,ex,
     &         slpd,pd,ssnmp)
               if(debug)write(*,*)'Called do_retrieval.'
            endif

c            if(ifm.eq.1) then   ! FOV center ray only
c            elseif(ifm.eq.3) then  ! Full numerical FOV integration (slow)
c            endif

c  Set error bars very small for "luft", otherwise it will
c  dominate the CO2/Luft ratio uncertainties.
            if(index(winfo,' luft ').gt.0 .and. nfp.gt.0) ex(1)=1.0e-08

c  Write .spt file for the first MSPT spectral fits
            if(ipfs.gt.0) then
               xfs=cx(ipfs)
            else
               xfs=0.0
            endif

            if(ipsg.gt.0) then
               xsg=cx(ipsg)
            else
               xsg=0.0
            endif

            if(freq_flag.eq.0) then
c               dopp=xfs*gint/((nus+nue)/2)  ! Telluric Reference Frame
               dopp=1.E-6*xfs  ! Telluric Reference Frame
            else
               dopp=1.E-6*(xsg+xfs+osds)  ! Solar reference frame
            endif
            apz=asza+zenoff
            if(ispec .le. mspt) then
               sptpath=sptfile(:lf-1)//specname_rl
               call write_spt(lunw_spt,flagso,sptpath,obsrvd,calcul,
     &         cont,cx,ex,startm+gint*xfs,dopp,effres,gint,overcol,pars,
     &         apz,obalt,zmin,xzo,rmsocl,peff,frac,pd,ssnmp,nmp,ntg,nfp)
            endif

cc  Solar_Gas_shift is the difference of the solar shift and the gas shift
cc  expressed in terms of the observed spectral point spacing
cc  Multiply by GINT to convert to cm-1.
cc   Divide by frequency to normalize to a dimensionless stretch which
cc  should then  be the same for all windows. Multiply by 10^6 for ppm.
c            if(index(winfo,' so ').gt.0) then
c               sgshift=1E6*2*gint/frqcen*
c     &         solar_gas_shift(cont_level,cont_tilt,obsrvd,calcul,ssnmp,nmp)
c            else
c               sgshift=0.0
c            endif

c  Write CBF and channel fringe coefficients to file. Make sure the
c  column for the spectrum name has enough room for the full spectrum
c  name plus two for quotes (to facilitate free form reading).
            cbffile_format = ''// 
     &       '(a40,1x,f8.5,f9.3,f7.3,1pe12.4,0p30f9.4)'
            write(cbffile_format(3:4),'(i2.2)') lspmax+2
            write(lunw_cbf,cbffile_format)
     &      '"'//specname_rl(:lspmax)//'"',
     &      wlimit(cfamp/cont_level, 'f8.5'),cfperiod*graw,cfphase,
     &      (cx(j+ntg),j=1,ncbf)

c  Write the state vector and uncertainties to the .col file
            call write_col(lunw_col,colfile_format,nspectra0-nspectra,
     &       specname_rl,lspmax,nit,xzo,xfs,cont_level,cont_tilt,
     &       cont_curv,rmsocl,xsg,zmin,oloscol,overcol,cx,ex,ntg,nfp)
c  And to the screen
            call write_col(6,colfile_format,nspectra0-nspectra,
     &       specname_rl,lspmax,nit,
     &       xzo,xfs,cont_level,cont_tilt,cont_curv,
     &       rmsocl,xsg,zmin,oloscol,overcol,cx,ex,ntg,nfp)

c            if(debug .and. nit.eq.mit+1) stop 'nit=mit+1' 

c  Output PD's (weighting functions), for subsequent use
c  in deriving averaging kernels.
            if( index(winfo,' ak ') .gt. 0) then
               open(lunw_ak,file=akpath(:lnbc(akpath))//'_'//
     &         specname_rl,status='unknown')
               if(debug)write(*,*)'Calling fm.....spxv(1,0)=',spxv(1)
               call fm(lunw_ak,slit,nhw,
     &         ifcsp,ifmsp,iptg,ipcl,ipfs,ipsg,ipzo,ipcf,
     &         ldec,spts,spxv,
     &         vac,splos,nlev_mav,ncp,rdec,
     &         cont_level,cont_tilt,cont_curv,xzo,
     &         cx,ntg,ncbf,nfp,cont,
     &         calcul,slpd,pd,nmp)
c               ynoise=2.5*cont_level*corrld/sngl(0.1d0+snr)
               ynoise=cont_level*corrld/sngl(0.1d0+snr) ! GCT 2020-01-11
c  Skip levels representing the cells, so start ilev at ncell+1.
               write(lunw_ak,*)(splos(ilev)*cp(ilev),ilev=ncell+1,
     &         nlev_ray)
               write(lunw_ak,*) pout/1013.25
               write(lunw_ak,*)(z(ilev),ilev=ncell+1,nlev_ray)
               write(lunw_ak,*)(p(ilev),ilev=ncell+1,nlev_ray)
               write(lunw_ak,*)(ynoise/apru(ifp),ifp=1,nfp)
               write(lunw_ak,*)((aprx(ifp)-cx(ifp))*ynoise/apru(ifp),
     &         ifp=1,nfp)
               close(lunw_ak)
            endif

c            if(rmsocl.le.0.0) rmsocl=0.000001
c  Use wlimit(rmsocl) to be consistent with collate_results
            rmspc=wlimit(dble(abs(100*rmsocl)),'f6.4')  ! rmsoclpc typically ~0.33%
c            wt= 1./(3*rmspc+1.0/3/rmspc)
            wt= 1./(3*rmspc+0.1)
            twt=twt+wt
            trmspc=trmspc+wt*rmspc
c            trmspc=trmspc+1.d0/
c     &      (wlimit(dble(abs(100*rmsocl)),'f6.4')/100+.0001)
            totit=totit+nit
            if(nfp.gt.0) then
               do jt=1,ntg
                  toterr(jt)=toterr(jt)+1.0/ex(jt)**2
                  tottc(jt)=tottc(jt)+cx(jt)/ex(jt)**2
                  tottc2(jt)=tottc2(jt)+((cx(jt)-1)/ex(jt))**2
               end do
            endif

         endif        !  (ncall.ge.1) then
         nspectra=nspectra+1
      end do   !  ispec=1,mspectra     Main fitting loop over spectra

      if(ncall.eq.0) then  ! Pre-secreening
         write(*,*) '# of spectra to be fitted = ',nspectra
         nspectra0=nspectra
      else  ! after fitting
         write(*,*)' Grid=',grid,'cm-1' 
         write(*,*)' Used ',float(nva+ncp)/mva,' of allocated memory'
         write(*,*)' Total number of iterations =',totit
         write(*,*)' Total number of spectra fitted  =',nspectra
         write(*,*)' Total number of mavblocks found =',mav_count
c         write(*,*)' Average % RMS fit =',100*avgrms/avgcl
         write(*,'(a,f8.5)')' Average % RMS fit =',trmspc/twt
         write(*,*)
         write(*,'(13x,9(4x,a,i1))') ('Target Gas ',jt,jt=1,min(ntg,9))
         write(*,'(a,9(f8.4,a2,f6.4))')' Average VSF =',
     &   (tottc(jt)/toterr(jt),'+-',
     &   sqrt(abs(tottc2(jt)/toterr(jt)-(tottc(jt)/toterr(jt)-1.)**2)),
     &   jt=1,min(ntg,9))
      endif
      close(lun_rlg)
      close(lunr_ray)
      close(lunr_mav)
      ncall=ncall+1
      return
      end


      function riair(w,t,p,h)
c  Edlen's formula for the refractive index of air at wavenumber w cm-1
c  T is the air temperature in degrees Celsius
c  P is the air pressure in mbar
C  H is the relative humidity in %
      real*8 w,t,p,h,pp,delt,hh,riair
c      F1F(w)=0.378125+w*w*(2.1414E-11+w*w*1.793E-21)
c      F2F(w)=0.0624-w*w*6.8E-14
      if(w.lt.0.0) w=0.0
      PP=(p-.3175+5.E-4*(p-745.)+.13*(t-10.))*.7500646
      DELT=PP*(1.-PP*(1.57E-8*t-1.049E-6))/(1.+3.661E-3*t)
      HH=h*EXP(1.52334+t*(.07217-2.9549E-4*t))/(100.+.3661*t)
      riair=1.D0+1.E-6*(DELT*(0.378125+w*w*(2.1414E-11+w*w*1.793E-21))
     & -HH*(0.0624-w*w*6.8E-14))
      return
      END


C***************************************************************************
C---------------------------------------------------------------------------
C***************************************************************************
      Subroutine CompAbs(sgmin,sgmax,dsig,T,Ptot,xCH4,Abs_HT)
C Subroutine computes the absorption coeffcients of CH4 diluted in air 
C in the 2NU3 R(6) manifold region
C
      Implicit None
      Integer nSigmx,nLinemx
      Integer Nsig,iexp
      Double precision Pi,A,aMass,T0,CC,sigC
      Double precision sgmin,sgmax,dsig,T,Ptot,xCH4
      Double precision Qr
c      Double precision Cte,Cte1,Cte2
      Double precision gi,QT0,QT
C Max number of Spectral points
      Parameter (nSigmx=800000)
c      Parameter (Nhitmx=2000)
      Parameter (Nlinemx=6)
C Hitran data
c      Double precision Sigmar,Srt0,S
c      Double precision HWairt0,HWselft0
c      Double precision Srt(Nhitmx),Shift(Nhitmx),Sigr(Nhitmx)
c      Double precision Alart(Nhitmx),GamD(Nhitmx)
c      Character Syml*5,vup*11,vlow*11
C 2nu3 R(6) manifold data
      integer iso,iline,Nline
      Double precision sig,strR6T0,EngR6,gam0,gamA,ngam0
      Double precision shift0,nshift0,Y,nY,nuVC,nnuVC
      Double precision gam2,ngam2,shift2,nshift2,eta
      Double precision sigR6(Nlinemx),strR6T(Nlinemx)
      Double precision gam0R6(Nlinemx),shift0R6(Nlinemx)
      Double precision Gamdop(Nlinemx),YR6(Nlinemx)
      Double precision nuVCR6(Nlinemx),gam2R6(Nlinemx)
      Double precision shift2R6(Nlinemx),etaR6(Nlinemx)

C Intermediate results
      Double precision Aloh(nSigmx),AbsR6(nSigmx),AbsR61(nSigmx)
      Double precision LS_pCqSDHC_R,LS_pCqSDHC_I
C Results (Absorption Coefficients)
      Double precision Abs_HT(nSigmx)
c      Double precision XX,YY,WR,WI
c      Common/CabsLM/Abs_HT(nSigmx)
C
C Constantes
C
      Pi=4.d0*datan(1.d0)
      A=1.4388d0
      aMass=16.0425d0
      T0=296.d0
      CC=0.1013/(1.38D-23*T)

      Nsig=int((sgmax-sgmin)/dsig)+1
      if (Nsig.gt.nSigmx) then
         write(6,*)'Set nSigmx at least equal to',Nsig
         stop
      endif
C
C*****Reading of HITRAN 2012 data base for the considered region********
C
c      Open(60,File='/home/toon/hatran/06_hit12.par',status='old')
c      ihit=0
cC Hit2012      
c68      Read(60,67,end=69)Iso,Sigmar,Srt0,HWairt0,HWselft0,Enr,
c     1Bt,S,vup,vlow,Ju,Jl,Syml
c67    Format(2x,i1,f12.6,e10.3,10x,f5.4,f5.3,f10.4,
c     1f4.2,f8.6,4x,a11,4x,a11,3x,i2,13x,i2,a5)    
c      if (Sigmar.lt.(sgmin-1.d0).or.sigmar.gt.(sgmax+1.d0)) go to 68
cC***************************************      
c      If (Iso.eq.1.and.vup.eq.'0 0 2 0  F2'.and.vlow.eq.'0 0 0 0 1A1'.
c     1and.Jl.eq.6) go to 68
c      ihit=ihit+1
cC Check that dimensions are fine. Stop if not
c      If ( ihit .GT. Nhitmx ) Then
c      Write(*,2000)                 
c      Stop
c      End If
c2000  Format(//,1x,'************ PROBLEM !!!! ******************',
c     /       /,1x,'Arrays in for Line data storage are too small',
c     /       /,1x,'raise the value of Nhitmx')
cC
c      Alart(ihit)=((T0/T)**Bt)*HWairt0
c      call BD_TIPS_2011(6,T0,Iso,gi,QT0)
c      call BD_TIPS_2011(6,T,Iso,gi,QT)
c      Qr=QT0/QT
c
c      Srt(ihit)=CC*Srt0*Qr*dexp(-A*Enr*(1./T-1./T0))
c     1  /(1.d0-dexp(-A*Sigmar/T0))/Sigmar
c      Sigr(ihit)=Sigmar
c      Shift(ihit)=S
c      GamD(ihit)=3.581163139d-7*DSQRT(T/aMass)*Sigr(ihit)
c      Goto 68                                  
c69    Close (60)
c      Nhit=ihit
cc      Write(6,*)'Number of small lines ',Nhit
cC**********************************************************************
c      do ihit=1,Nhit
c      Sigr(ihit)=Sigr(ihit)+Shift(ihit)*Ptot
c      end do
c      Do 170 iexp=1,nSig
c      Aloh(iexp)=0.d0
c170    Continue
c
c      do 63 ihit=1,Nhit
c      Do 62 iexp=1,Nsig
c      SigC=sgmin+(iexp-1)*DSig
c      Cte2=SigC*(1.d0-DEXP(-1.4388d0*SigC/T))
c      
c      Cte=DSQRT(DLOG(2.d0))/GamD(ihit)
c      Cte1=Cte/DSQRT(Pi)
c
c      XX=(Sigr(ihit)-SigC)*Cte
c      YY=Alart(ihit)*Ptot*Cte
c      Call CPF(XX,YY,WR,WI)
c
c      aa=WR*Cte1*Srt(ihit)
c      Aloh(iexp)=Aloh(iexp)+aa*Cte2
c62      Continue
c63      Continue
cC
C*******Reading data of the 2nu3 R(6) manifold********************
      open(9,file='/home/toon/hatran/Parameters_fit_det.par')      
      read(9,*)
      read(9,*)
      iline=0
90    read(9,901,end=19)iso,sig,strR6T0,EngR6,gam0,gamA,ngam0,
     1shift0,nshift0,Y,nY,nuVC,nnuVC,gam2,ngam2,shift2,nshift2,eta
901   format(2x,i1,f12.6,e12.5,f10.4,e15.8,f8.5,e14.6,2x,
     1e15.8,12x,f12.8,e17.8,e15.6,3(e17.8,f10.5),e17.8)

      iline=iline+1
      call BD_TIPS_2011(6,T0,Iso,gi,QT0)
      call BD_TIPS_2011(6,T,Iso,gi,QT)
      Qr=QT0/QT
      StrR6T(iline)=CC*strR6T0*Qr*dexp(-A*EngR6*(1./T-1./T0))*
     1(1.d0-dexp(-A*Sig/T))/(1.d0-dexp(-A*Sig/T0))
      SigR6(iline)=Sig
      gam0R6(iline)=Ptot*gam0*(T0/T)**ngam0
c     shift0R6(iline)=Ptot*shift0*(T0/T)**nshift0         
      shift0R6(iline)=Ptot*(shift0+nshift0*(T-T0))   
        
      YR6(iline)=Ptot*Y*(T0/T)**nY
      nuVCR6(iline)=Ptot*nuVC*(T0/T)**nnuVC
      gam2R6(iline)=Ptot*gam2*(T0/T)**ngam2
      shift2R6(iline)=Ptot*shift2*(T0/T)**nshift2
      etaR6(iline)=eta
      GamDop(iline)=3.581163139d-7*DSQRT(T/aMass)*Sig
      go to 90
19    close(9)
      Nline=iline

C**********************************************************************
      do iexp=1,nSig
         AbsR6(iexp)=0.d0
         AbsR61(iexp)=0.d0
      end do

      do iline=1,Nline
         do iexp=1,Nsig
            SigC=sgmin+(iexp-1)*DSig
            call pCqSDHC(sigR6(iline),GamDop(iline),gam0R6(iline),
     1      Gam2R6(iline),Shift0R6(iline),Shift2R6(iline),nuVCR6(iline),
     2      etaR6(iline),SigC,LS_pCqSDHC_R,LS_pCqSDHC_I)
      
C -Y for the new parameters 
            AbsR6(iexp)=AbsR6(iexp)+
     &     StrR6T(iline)*(LS_pCqSDHC_R-YR6(iline)*LS_pCqSDHC_I)
         end do
      end do  ! iline=1,Nline

      Do iexp=1,Nsig
         Abs_HT(iexp)=(AbsR6(iexp)+Aloh(iexp))*Ptot*xCH4
      end do

      Return
      End Subroutine CompAbs

!****************************************************************************      
C******************************************************
C***************************************************

      subroutine pCqSDHC(sg0,GamD,Gam0,Gam2,Shift0,Shift2,anuVC,eta,
     &sg,LS_pCqSDHC_R,LS_pCqSDHC_I)
C-------------------------------------------------
C      "pCqSDHC": partially-Correlated quadratic-Speed-Dependent Hard-Collision
C      Subroutine to Compute the complex normalized spectral shape of an 
C      isolated line by the pCqSDHC model following the two references:
C      [1] Ngo NH, Lisak D, Tran H, Hartmann J-M. An isolated line-shape model
C      to go beyond the Voigt profile in spectroscopic databases and radiative 
C      transfer codes. J Quant Radiat Transfer 2013;129:89-100.
C      [2] Tran H, Ngo NH, Hartmann J-M. Efficient computation of some speed-dependent 
C      isolated line profiles. J Quant Radiat Transfer 2013;129:199-203.
C
C      Input/Output Parameters of Routine (Arguments or Common)
C      ---------------------------------
C      T:      Temperature in Kelvin (Input).
C      amM1:   Molar mass of the absorber in g/mol(Input).
C      sg0:    Unperturbed line position in cm-1 (Input).
C      GamD:   Doppler HWHM in cm-1 (Input)
C      Gam0:   Speed-averaged line-width in cm-1 (Input). 
C      Gam2:   Speed dependence of the line-width in cm-1 (Input).
C      anuVC:  Velocity-changing frequency in cm-1 (Input).
C      eta:    Correlation parameter, No unit (Input).
C      Shift0: Speed-averaged line-shift in cm-1 (Input).
C      Shift2: Speed dependence of the line-shift in cm-1 (Input) 
C      sg:     Current WaveNumber of the Computation in cm-1 (Input).
C
C      Output Quantities (through Common Statements)
C      -----------------
C      LS_pCqSDHC_R: Real part of the normalized spectral shape (cm)
C      LS_pCqSDHC_I: Imaginary part of the normalized spectral shape (cm)
C
C      Called Routines: 'CPF' (Complex Probability Function)
C      ---------------  'CPF3' (Complex Probability Function for the region 3)
C
C      Called By: Main Program
C      ---------
C
C      Double Precision Version
C
C-------------------------------------------------
      implicit none
      double precision sg0,GamD
      double precision Gam0,Gam2,anuVC,eta,Shift0,Shift2
      double precision sg
      double precision pi,rpi,cte
      double precision xz1,xz2,yz1,yz2,xXb,yXb
      double precision wr1,wi1,wr2,wi2,wrb,wib
      double precision SZ1,SZ2,DSZ,SZmx,SZmn
      double precision LS_pCqSDHC_R,LS_pCqSDHC_I
      double complex c0,c2,c0t,c2t
      double complex X,Y,iz,Z1,Z2,csqrtY
      double complex Aterm,Bterm,LS_pCqSDHC
C
C-------------------------------------------------
C
      cte=dsqrt(dlog(2.D0))/GamD
      pi=4.d0*datan(1.d0)
      rpi=dsqrt(pi)
      iz=dcmplx(0.d0,1.d0)
c Calculating the different parameters 
      c0=dcmplx(Gam0,Shift0)
      c2=dcmplx(Gam2,Shift2)
      c0t=(1.d0-eta)*(c0-1.5d0*c2)+anuVC
      c2t=(1.d0-eta)*c2
C
c      if (cdabs(c2t).eq.0.d0) go to 110 ! GCT 2020-01-26
      if (cdabs(c2t).le.0.d0) go to 110
      X=(iz*(sg0-sg)+c0t)/c2t
      Y=1.d0/((2.d0*cte*c2t))**2
      csqrtY=(Gam2-iz*Shift2)/(2.d0*cte*(1.d0-eta)*(Gam2**2+Shift2**2))
      if (cdabs(X).le.3.d-8*cdabs(Y)) go to 120
      if (cdabs(Y).le.1.d-15*cdabs(X)) go to 140
c calculating Z1 and Z2
      Z1=cdsqrt(X+Y)-csqrtY
      Z2=Z1+2.d0*csqrtY
c calculating the real and imaginary parts of Z1 and Z2
      xZ1=-dimag(Z1)
      yZ1=dreal(Z1)
      xZ2=-dimag(Z2)
      yZ2=dreal(Z2)
c check if Z1 and Z2 are close to each other
      SZ1=dsqrt(xZ1*xZ1+yZ1*yZ1)
      SZ2=dsqrt(xZ2*xZ2+yZ2*yZ2)
      DSZ=dabs(SZ1-SZ2)
      SZmx=dmax1(SZ1,SZ2)
      SZmn=dmin1(SZ1,SZ2)
c when Z1 and Z2 are close to each other, ensure that they are in 
c the same interval of CPF 
      if (DSZ.le.1.d0.and.SZmx.gt.8.d0.and.SZmn.le.8.d0) then
         Call CPF3 ( xZ1, yZ1, WR1, WI1 ) 
         Call CPF3 ( xZ2, yZ2, WR2, WI2 ) 
      else
         Call CPF ( xZ1, yZ1, WR1, WI1 ) 
         Call CPF ( xZ2, yZ2, WR2, WI2 ) 
      endif
c calculating the A and B terms of the profile
      Aterm=rpi*cte*(dcmplx(wr1,wi1)-dcmplx(wr2,wi2))
      Bterm=(-1.d0+
     ,rpi/(2.d0*csqrtY)*(1.d0-Z1**2)*dcmplx(wr1,wi1)-
     ,rpi/(2.d0*csqrtY)*(1.d0-Z2**2)*dcmplx(wr2,wi2))/
     ,c2t
      go to 10
c when c2t=0
110   continue
      Z1=(iz*(sg0-sg)+c0t)*cte
      xZ1=-dimag(Z1)
      yZ1=dreal(Z1)
      Call CPF ( xZ1, yZ1, WR1, WI1 )
      Aterm=rpi*cte*dcmplx(WR1,WI1)
      if (cdabs(Z1).le.4.d3) then
         Bterm=rpi*cte*((1.d0-Z1**2)*dcmplx(WR1,WI1)+Z1/rpi)
      else
         Bterm=cte*
     &   (rpi*dcmplx(WR1,WI1)+0.5d0/Z1-0.75d0/(Z1**3))
      endif
      go to 10
c when abs(Y) is much larger than abs(X)
120   continue
      Z1=(iz*(sg0-sg)+c0t)*cte
      Z2=cdsqrt(X+Y)+csqrtY
      xZ1=-dimag(z1)
      yZ1=dreal(z1)
      xZ2=-dimag(z2)
      yZ2=dreal(z2)
      Call CPF ( xZ1, yZ1, WR1, WI1 )
      Call CPF ( xZ2, yZ2, WR2, WI2 ) 
      Aterm=rpi*cte*(dcmplx(WR1,WI1)-dcmplx(WR2,WI2))
      Bterm=(-1.d0+
     ,rpi/(2.d0*csqrtY)*(1.d0-Z1**2)*dcmplx(wr1,wi1)-
     ,rpi/(2.d0*csqrtY)*(1.d0-Z2**2)*dcmplx(wr2,wi2))/
     ,c2t
      go to 10
c when abs(X) is much larger than abs(Y)
140   continue
      xZ1=-dimag(cdsqrt(X+Y))
      yZ1=dreal(cdsqrt(X+Y))
      Call CPF ( xZ1, yZ1, WR1, WI1 ) 
      if (cdabs(cdsqrt(X)).le.4.d3) then
         xXb=-dimag(cdsqrt(X))
         yXb=dreal(cdsqrt(X))
         Call CPF ( xXb, yXb, WRb, WIb ) 
         Aterm=(2.d0*rpi/c2t)*(1.d0/rpi-cdsqrt(X)*dcmplx(WRb,WIb))
         Bterm=(1.d0/c2t)*(-1.d0+
     ,   2.d0*rpi*(1.d0-X-2.d0*Y)*(1.d0/rpi-cdsqrt(X)*dcmplx(wrb,wib))+
     ,   2.d0*rpi*cdsqrt(X+Y)*dcmplx(wr1,wi1))
cc and when abs(X) is much larger than 1
      else
         Aterm=(1.d0/c2t)*(1.d0/X-1.5d0/(X**2))
         Bterm=(1.d0/c2t)*(-1.d0+(1.d0-X-2.d0*Y)*
     ,   (1.d0/X-1.5d0/(X**2))+
     ,   2.d0*rpi*cdsqrt(X+Y)*dcmplx(wr1,wi1))
      endif
c
10    continue
c
      LS_pCqSDHC=(1.d0/pi)*(Aterm/
     ,(1.d0-(anuVC-eta*(C0-1.5d0*C2))*Aterm+
     ,eta*C2*Bterm))

      LS_pCqSDHC_R=dreal(LS_pCqSDHC)
      LS_pCqSDHC_I=dimag(LS_pCqSDHC)

   
      Return
      End Subroutine pCqSDHC
C
C-------------------------------------------------
C

      Subroutine CPF(X,Y,WR,WI)
C-------------------------------------------------
C "CPF": Complex Probability Function
C .........................................................
C         .       Subroutine to Compute the Complex       .
C         .        Probability Function W(z=X+iY)         .
C         .     W(z)=exp(-z**2)*Erfc(-i*z) with Y>=0      .
C         .    Which Appears when Convoluting a Complex   .
C         .     Lorentzian Profile by a Gaussian Shape    .
C         .................................................
C
C             WR : Real Part of W(z)
C             WI : Imaginary Part of W(z)
C
C This Routine was Taken from the Paper by J. Humlicek, which 
C is Available in Page 309 of Volume 21 of the 1979 Issue of
C the Journal of Quantitative Spectroscopy and Radiative Transfer
C Please Refer to this Paper for More Information
C
C Accessed Files:  None
C --------------
C
C Called Routines: None                               
C ---------------                                 
C
C Called By: 'CompAbs' (COMPute ABSorpton)
C ---------
C
C Double Precision Version
C
C-------------------------------------------------
C      
      Implicit None
      Integer I
      double complex zm1,zm2,zterm,zsum,zone,zi
      Double Precision X,Y,WR,WI
      Double Precision T,U,S,Y1,Y2,Y3,R,R2,D,D1,D2,D3,D4
      Double Precision TT(15),pipwoeronehalf
C      
      Dimension T(6),U(6),S(6)
      Data T/.314240376d0,.947788391d0,1.59768264d0,2.27950708d0,
     & 3.02063703d0,3.8897249d0/
      Data U/1.01172805d0,-.75197147d0,1.2557727d-2,1.00220082d-2,
     & -2.42068135d-4,5.00848061d-7/
      Data S/1.393237d0,.231152406d0,-.155351466d0,6.21836624d-3,
     & 9.19082986d-5,-6.27525958d-7/
      Data zone,zi/(1.d0,0.D0),(0.d0,1.D0)/
      data tt/0.5d0,1.5d0,2.5d0,3.5d0,4.5d0,5.5d0,6.5d0,7.5d0,8.5d0,
     & 9.5d0,10.5d0,11.5d0,12.5d0,13.5d0,14.5d0/
      data pipwoeronehalf/0.564189583547756d0/

C new Region 3
      if(dsqrt(x*x+y*Y).gt.8.D0)then
         zm1=zone/dcmplx(x,y)
         zm2=zm1*zm1
         zsum=zone
         zterm=zone
         do i=1,15
            zterm=zterm*zm2*tt(i)
            zsum=zsum+zterm
         end do
         zsum=zsum*zi*zm1*pipwoeronehalf
         wr=dreal(zsum)
         wi=dimag(zsum)
         return
      end if
C
      WR=0.d0
      WI=0.d0
      Y1=Y+1.5d0
      Y2=Y1*Y1
      If( (Y.GT.0.85d0) .OR. (DABS(X).LT.(18.1d0*Y+1.65d0)) )GoTo 2
C
C       Region 2
C
      If( DABS(X).LT.12.d0 )WR=DEXP(-X*X)
      Y3=Y+3.d0
      Do I=1,6
         R=X-T(I)
         R2=R*R
         D=1.d0/(R2+Y2)
         D1=Y1*D
         D2=R*D
         WR=WR+Y*(U(I)*(R*D2-1.5d0*D1)+S(I)*Y3*D2)/(R2+2.25d0)
         R=X+T(I)
         R2=R*R
         D=1.d0/(R2+Y2)
         D3=Y1*D
         D4=R*D
         WR=WR+Y*(U(I)*(R*D4-1.5d0*D3)-S(I)*Y3*D4)/(R2+2.25d0)
         WI=WI+U(I)*(D2+D4)+S(I)*(D1-D3)
      end do
      Return
C
C    Region 1
C
 2    Continue
      Do I=1,6
         R=X-T(I)
         D=1.d0/(R*R+Y2)
         D1=Y1*D
         D2=R*D
         R=X+T(I)
         D=1.d0/(R*R+Y2)
         D3=Y1*D
         D4=R*D
         WR=WR+U(I)*(D1+D3)-S(I)*(D2-D4)
         WI=WI+U(I)*(D2+D4)+S(I)*(D1-D3)
      end do
      Return
      End

      Subroutine CPF3(X,Y,WR,WI)
C-------------------------------------------------
C "CPF": Complex Probability Function
C .........................................................
C         .       Subroutine to Compute the Complex       .
C         .        Probability Function W(z=X+iY)         .
C         .     W(z)=exp(-z**2)*Erfc(-i*z) with Y>=0      .
C         .    Which Appears when Convoluting a Complex   .
C         .     Lorentzian Profile by a Gaussian Shape    .
C         .................................................
C
C             WR : Real Part of W(z)
C             WI : Imaginary Part of W(z)
C
C This Routine takes into account the region 3 only, i.e. when sqrt(x**2+y**2)>8. 
C
C Accessed Files:  None
C --------------
C
C Called Routines: None                               
C ---------------                                 
C
C Called By: 'pCqSDHC'
C ---------
C
C Double Precision Version
C 
C-------------------------------------------------
C      
      Implicit None
      Integer I
      double complex zm1,zm2,zterm,zsum,zone,zi
      Double Precision X,Y,WR,WI
      Double Precision TT(15),pipwoeronehalf
C      
      Data zone,zi/(1.d0,0.D0),(0.d0,1.D0)/
      data tt/0.5d0,1.5d0,2.5d0,3.5d0,4.5d0,5.5d0,6.5d0,7.5d0,8.5d0,
     &        9.5d0,10.5d0,11.5d0,12.5d0,13.5d0,14.5d0/
      data pipwoeronehalf/0.564189583547756d0/

C Region 3
      zm1=zone/dcmplx(x,y)
      zm2=zm1*zm1
      zsum=zone
      zterm=zone
      do i=1,15
         zterm=zterm*zm2*tt(i)
         zsum=zsum+zterm
      end do
      zsum=zsum*zi*zm1*pipwoeronehalf
      wr=dreal(zsum)
      wi=dimag(zsum)
      return
      End
C      --

         
C***********************
      Subroutine BD_TIPS_2011(
     I MOL,     ! HITRAN molecule number
     I Temp,    ! temperature in K
     I ISO,     ! isotopomer index
     O gi,      ! state independent degeneracy factor
     O QT)      ! total internal partition sum
C*********************** 
C    Subroutine BD_TIPS_2010 written by R.R. Gamache
C 
c
c...date last changed 17         June, 2011
c  --  UPDATES  --
C    18 March 2011 addition of molecules 43-51, 838,837 of CO2, 
C     1222 of C2H2, CH3D, 13CH3D, 13C12CH6,CH3Br
C    19 July 2010 correction of the sign of D in CO constants
C    13 July 2010 addition of several new species/isotopologues: ^13C^18O2, 
C       ^18O^13C^17O, ^13CH3D,  ^13C12CH2, CF4: note a number of the block 
C       data codes were updated as well
C    15 May 2008 Dijon Q(T) values for CH4
C    18 December 2003 better vibrational fundamentals for PH3
c
C    This program calculates the total internal
c    partition sum (TIPS) for a given molecule, isotopic species, and
c    temperature.  Current limitations are the molecular species on the
c    HITRAN molecular database and the temperature range 70 - 3000 K.
c
c...This program calculates the TIPS by 4-point LaaGrange interpolation
c
c..  JQSRT - 82, 401-412, 2003
c..  J. Fischer(a) R.R. Gamache(a&), A. Goldman(b),L.S. Rothman(c), and A. Perrin(d)
c..  
c..  (a)  Department of Environmental, Earth, and Atmospheric Sciences, 
c..       University of Massachusetts Lowell, Lowell, MA 01854, U.S.A.
c..  
c..  (b)  Department of Physics, University of Denver, Denver, CO 80208, U.S.A.
c..  
c..  (c)  Harvard-Smithsonian Center for Astrophysics, 60 Garden St, Cambridge, MA 02138 USA
c..  
c..  (d)  Laboratoire de Photophysique Moléculaire, Université Paris Sud, 91405 Orsay, FRANCE
c..  
c..  &  Corresponding author. Email address: Robert_Gamache@uml.edu
c..  Abstract
c..        Total internal partition sums (TIPS) are calculated for all molecular species in 
c..  the 2000 HITRAN database.  In addition, the TIPS for 13 other isotopomers/isotopologues 
c..  of ozone and carbon dioxide are presented.  The calculations address the corrections 
c..  suggested by Goldman et al. (JQSRT 2000;66:55-86).  The calculations consider the 
c..  temperature range 70-3000 K to be applicable to a variety of remote sensing needs.  
c..  The method of calculation for each molecular species is stated and comparisons with 
c..  data from the literature are discussed.  A new method of recall for the partition sums, 
c..  Lagrange 4-point interpolation, is developed.  This method, unlike previous versions of 
c..  the TIPS code, allows all molecular species to be considered.  
c
      implicit DOUBLE PRECISION (a-h,o-z)
 
c++
c      INCLUDE 'Species_2011.cmn'
c++
      PARAMETER (NT=119)
c++
c      INCLUDE 'Isotops.cmn'
c      INCLUDE 'MOLEC.cmn'
c++:  bd-QT
      COMMON/Temperatures/tdat(NT)
 
      data Tdat/  60.,  85., 110., 135., 160., 185., 210., 235.,
     + 260., 285., 310., 335., 360., 385., 410., 435., 460., 485.,
     + 510., 535., 560., 585., 610., 635., 660., 685., 710., 735.,
     + 760., 785., 810., 835., 860., 885., 910., 935., 960., 985.,
     +1010.,1035.,1060.,1085.,1110.,1135.,1160.,1185.,1210.,1235.,
     +1260.,1285.,1310.,1335.,1360.,1385.,1410.,1435.,1460.,1485.,
     +1510.,1535.,1560.,1585.,1610.,1635.,1660.,1685.,1710.,1735.,
     +1760.,1785.,1810.,1835.,1860.,1885.,1910.,1935.,1960.,1985.,
     +2010.,2035.,2060.,2085.,2110.,2135.,2160.,2185.,2210.,2235.,
     +2260.,2285.,2310.,2335.,2360.,2385.,2410.,2435.,2460.,2485.,
     +2510.,2535.,2560.,2585.,2610.,2635.,2660.,2685.,2710.,2735.,
     +2760.,2785.,2810.,2835.,2860.,2885.,2910.,2935.,2960.,2985.,
     +3010./
C
C     GO TO PARTICULAR MOLECULE:  
 
      IF(MOL.EQ.1) THEN 
         CALL QT_H2O(Temp,ISO,gi,QT)
         go to 100
      ENDIF 
C 
      IF(MOL.EQ.2) THEN 
         CALL QT_CO2(Temp,ISO,gi,QT)
         go to 100
      ENDIF 
C
      IF(MOL.EQ.6) THEN 
         CALL QT_CH4(Temp,ISO,gi,QT)
         go to 100
      ENDIF 

 100  RETURN
      END 
c
c     *****************
      Subroutine QT_H2O   (                       
     I T,       ! temperature in K 
     I iso,       ! isotope code (HITRAN INDEX)
     O gsi,       ! state independent nuclear degeneracyfactor
     O QT)       ! Total Internal Partition Function
 
      implicit DOUBLE PRECISION (a-h,o-z)
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 6), QofT( 6,119),Q(NT)
      data xgj/ 1.,1.,6.,6.,6.,36/
c...      H2O
c...        --       161
      data (QofT( 1,J),J=1,119)/ 0.16824E+02, 0.27771E+02, 0.40408E+02,
     + 0.54549E+02, 0.70054E+02, 0.86817E+02, 0.10475E+03, 0.12380E+03,
     + 0.14391E+03, 0.16503E+03, 0.18714E+03, 0.21021E+03, 0.23425E+03,
     + 0.25924E+03, 0.28518E+03, 0.31209E+03, 0.33997E+03, 0.36883E+03,
     + 0.39870E+03, 0.42959E+03, 0.46152E+03, 0.49452E+03, 0.52860E+03,
     + 0.56380E+03, 0.60015E+03, 0.63766E+03, 0.67637E+03, 0.71631E+03,
     + 0.75750E+03, 0.79999E+03, 0.84380E+03, 0.88897E+03, 0.93553E+03,
     + 0.98353E+03, 0.10330E+04, 0.10840E+04, 0.11365E+04, 0.11906E+04,
     + 0.12463E+04, 0.13037E+04, 0.13628E+04, 0.14237E+04, 0.14863E+04,
     + 0.15509E+04, 0.16173E+04, 0.16856E+04, 0.17559E+04, 0.18283E+04,
     + 0.19028E+04, 0.19793E+04, 0.20581E+04, 0.21391E+04, 0.22224E+04,
     + 0.23080E+04, 0.24067E+04, 0.24975E+04, 0.25908E+04, 0.26867E+04,
     + 0.27853E+04, 0.28865E+04, 0.29904E+04, 0.30972E+04, 0.32068E+04,
     + 0.33194E+04, 0.34349E+04, 0.35535E+04, 0.36752E+04, 0.38001E+04,
     + 0.39282E+04, 0.40597E+04, 0.41945E+04, 0.43327E+04, 0.44745E+04,
     + 0.46199E+04, 0.47688E+04, 0.49215E+04, 0.50780E+04, 0.52384E+04,
     + 0.54027E+04, 0.55710E+04, 0.57434E+04, 0.59200E+04, 0.61008E+04,
     + 0.62859E+04, 0.64754E+04, 0.66693E+04, 0.68679E+04, 0.70710E+04,
     + 0.72788E+04, 0.74915E+04, 0.77090E+04, 0.79315E+04, 0.81590E+04,
     + 0.83917E+04, 0.86296E+04, 0.88728E+04, 0.91214E+04, 0.93755E+04,
     + 0.96351E+04, 0.99005E+04, 0.10171E+05, 0.10448E+05, 0.10731E+05,
     + 0.11020E+05, 0.11315E+05, 0.11617E+05, 0.11924E+05, 0.12238E+05,
     + 0.12559E+05, 0.12886E+05, 0.13220E+05, 0.13561E+05, 0.13909E+05,
     + 0.14263E+05, 0.14625E+05, 0.14995E+05, 0.15371E+05, 0.15755E+05,
     + 0.16147E+05/
c...        --       181
      data (QofT( 2,J),J=1,119)/ 0.15960E+02, 0.26999E+02, 0.39743E+02,
     + 0.54003E+02, 0.69639E+02, 0.86543E+02, 0.10463E+03, 0.12384E+03,
     + 0.14412E+03, 0.16542E+03, 0.18773E+03, 0.21103E+03, 0.23531E+03,
     + 0.26057E+03, 0.28681E+03, 0.31406E+03, 0.34226E+03, 0.37130E+03,
     + 0.40135E+03, 0.43243E+03, 0.46456E+03, 0.49777E+03, 0.53206E+03,
     + 0.56748E+03, 0.60405E+03, 0.64179E+03, 0.68074E+03, 0.72093E+03,
     + 0.76238E+03, 0.80513E+03, 0.84922E+03, 0.89467E+03, 0.94152E+03,
     + 0.98982E+03, 0.10396E+04, 0.10909E+04, 0.11437E+04, 0.11982E+04,
     + 0.12543E+04, 0.13120E+04, 0.13715E+04, 0.14328E+04, 0.14959E+04,
     + 0.15608E+04, 0.16276E+04, 0.16964E+04, 0.17672E+04, 0.18401E+04,
     + 0.19151E+04, 0.19922E+04, 0.20715E+04, 0.21531E+04, 0.22370E+04,
     + 0.23232E+04, 0.24118E+04, 0.25030E+04, 0.25967E+04, 0.26929E+04,
     + 0.27918E+04, 0.28934E+04, 0.29978E+04, 0.31050E+04, 0.32151E+04,
     + 0.33281E+04, 0.34441E+04, 0.35632E+04, 0.36854E+04, 0.38108E+04,
     + 0.39395E+04, 0.40715E+04, 0.42070E+04, 0.43459E+04, 0.44883E+04,
     + 0.46343E+04, 0.47840E+04, 0.49374E+04, 0.50946E+04, 0.52558E+04,
     + 0.54209E+04, 0.55900E+04, 0.57632E+04, 0.59407E+04, 0.61224E+04,
     + 0.63084E+04, 0.64988E+04, 0.66938E+04, 0.68933E+04, 0.70975E+04,
     + 0.73064E+04, 0.75202E+04, 0.77389E+04, 0.79625E+04, 0.81913E+04,
     + 0.84252E+04, 0.86644E+04, 0.89089E+04, 0.91588E+04, 0.94143E+04,
     + 0.96754E+04, 0.99422E+04, 0.10215E+05, 0.10493E+05, 0.10778E+05,
     + 0.11068E+05, 0.11365E+05, 0.11668E+05, 0.11977E+05, 0.12293E+05,
     + 0.12616E+05, 0.12945E+05, 0.13281E+05, 0.13624E+05, 0.13973E+05,
     + 0.14330E+05, 0.14694E+05, 0.15066E+05, 0.15445E+05, 0.15831E+05,
     + 0.16225E+05/
c...        --       171
      data (QofT( 3,J),J=1,119)/ 0.95371E+02, 0.16134E+03, 0.23750E+03,
     + 0.32273E+03, 0.41617E+03, 0.51722E+03, 0.62540E+03, 0.74036E+03,
     + 0.86185E+03, 0.98970E+03, 0.11238E+04, 0.12642E+04, 0.14097E+04,
     + 0.15599E+04, 0.17159E+04, 0.18777E+04, 0.20453E+04, 0.22188E+04,
     + 0.23983E+04, 0.25840E+04, 0.27760E+04, 0.29743E+04, 0.31792E+04,
     + 0.33907E+04, 0.36091E+04, 0.38346E+04, 0.40672E+04, 0.43072E+04,
     + 0.45547E+04, 0.48100E+04, 0.50732E+04, 0.53446E+04, 0.56244E+04,
     + 0.59128E+04, 0.62100E+04, 0.65162E+04, 0.68317E+04, 0.71567E+04,
     + 0.74915E+04, 0.78363E+04, 0.81914E+04, 0.85571E+04, 0.89335E+04,
     + 0.93211E+04, 0.97200E+04, 0.10131E+05, 0.10553E+05, 0.10988E+05,
     + 0.11435E+05, 0.11895E+05, 0.12368E+05, 0.12855E+05, 0.13356E+05,
     + 0.13870E+05, 0.14399E+05, 0.14943E+05, 0.15502E+05, 0.16076E+05,
     + 0.16666E+05, 0.17272E+05, 0.17895E+05, 0.18534E+05, 0.19191E+05,
     + 0.19865E+05, 0.20557E+05, 0.21267E+05, 0.21996E+05, 0.22744E+05,
     + 0.23512E+05, 0.24299E+05, 0.25106E+05, 0.25935E+05, 0.26784E+05,
     + 0.27655E+05, 0.28547E+05, 0.29462E+05, 0.30400E+05, 0.31361E+05,
     + 0.32345E+05, 0.33353E+05, 0.34386E+05, 0.35444E+05, 0.36527E+05,
     + 0.37637E+05, 0.38772E+05, 0.39934E+05, 0.41124E+05, 0.42341E+05,
     + 0.43587E+05, 0.44861E+05, 0.46165E+05, 0.47498E+05, 0.48862E+05,
     + 0.50256E+05, 0.51682E+05, 0.53139E+05, 0.54629E+05, 0.56152E+05,
     + 0.57708E+05, 0.59299E+05, 0.60923E+05, 0.62583E+05, 0.64279E+05,
     + 0.66011E+05, 0.67779E+05, 0.69585E+05, 0.71429E+05, 0.73312E+05,
     + 0.75234E+05, 0.77195E+05, 0.79197E+05, 0.81240E+05, 0.83325E+05,
     + 0.85452E+05, 0.87622E+05, 0.89835E+05, 0.92093E+05, 0.94395E+05,
     + 0.96743E+05/
c...        --       162
      data (QofT( 4,J),J=1,119)/ 0.75792E+02, 0.12986E+03, 0.19244E+03,
     + 0.26253E+03, 0.33942E+03, 0.42259E+03, 0.51161E+03, 0.60619E+03,
     + 0.70609E+03, 0.81117E+03, 0.92132E+03, 0.10365E+04, 0.11567E+04,
     + 0.12820E+04, 0.14124E+04, 0.15481E+04, 0.16891E+04, 0.18355E+04,
     + 0.19876E+04, 0.21455E+04, 0.23092E+04, 0.24791E+04, 0.26551E+04,
     + 0.28376E+04, 0.30268E+04, 0.32258E+04, 0.34288E+04, 0.36392E+04,
     + 0.38571E+04, 0.40828E+04, 0.43165E+04, 0.45584E+04, 0.48089E+04,
     + 0.50681E+04, 0.53363E+04, 0.56139E+04, 0.59009E+04, 0.61979E+04,
     + 0.65049E+04, 0.68224E+04, 0.71506E+04, 0.74898E+04, 0.78403E+04,
     + 0.82024E+04, 0.85765E+04, 0.89628E+04, 0.93618E+04, 0.97736E+04,
     + 0.10199E+05, 0.10637E+05, 0.11090E+05, 0.11557E+05, 0.12039E+05,
     + 0.12535E+05, 0.13047E+05, 0.13575E+05, 0.14119E+05, 0.14679E+05,
     + 0.15257E+05, 0.15851E+05, 0.16464E+05, 0.17094E+05, 0.17743E+05,
     + 0.18411E+05, 0.19098E+05, 0.19805E+05, 0.20532E+05, 0.21280E+05,
     + 0.22049E+05, 0.22840E+05, 0.23652E+05, 0.24487E+05, 0.25345E+05,
     + 0.26227E+05, 0.27132E+05, 0.28062E+05, 0.29016E+05, 0.29997E+05,
     + 0.31002E+05, 0.32035E+05, 0.33094E+05, 0.34180E+05, 0.35295E+05,
     + 0.36438E+05, 0.37610E+05, 0.38812E+05, 0.40044E+05, 0.41306E+05,
     + 0.42600E+05, 0.43926E+05, 0.45284E+05, 0.46675E+05, 0.48100E+05,
     + 0.49559E+05, 0.51053E+05, 0.52583E+05, 0.54148E+05, 0.55750E+05,
     + 0.57390E+05, 0.59067E+05, 0.60783E+05, 0.62539E+05, 0.64334E+05,
     + 0.66170E+05, 0.68047E+05, 0.69967E+05, 0.71929E+05, 0.73934E+05,
     + 0.75983E+05, 0.78078E+05, 0.80217E+05, 0.82403E+05, 0.84636E+05,
     + 0.86917E+05, 0.89246E+05, 0.91625E+05, 0.94053E+05, 0.96533E+05,
     + 0.99064E+05/
c...        --       182
      data (QofT( 5,J),J=1,119)/ 0.82770E+02, 0.13749E+03, 0.20083E+03,
     + 0.27176E+03, 0.34955E+03, 0.43370E+03, 0.52376E+03, 0.61944E+03,
     + 0.72050E+03, 0.82679E+03, 0.93821E+03, 0.10547E+04, 0.11763E+04,
     + 0.13031E+04, 0.14350E+04, 0.15723E+04, 0.17150E+04, 0.18633E+04,
     + 0.20172E+04, 0.21770E+04, 0.23429E+04, 0.25149E+04, 0.26934E+04,
     + 0.28784E+04, 0.30702E+04, 0.32690E+04, 0.34750E+04, 0.36885E+04,
     + 0.39096E+04, 0.41386E+04, 0.43758E+04, 0.46213E+04, 0.48755E+04,
     + 0.51386E+04, 0.54109E+04, 0.56927E+04, 0.59841E+04, 0.62856E+04,
     + 0.65973E+04, 0.69197E+04, 0.72529E+04, 0.75973E+04, 0.79533E+04,
     + 0.83210E+04, 0.87009E+04, 0.90933E+04, 0.94985E+04, 0.99168E+04,
     + 0.10348E+05, 0.10794E+05, 0.11254E+05, 0.11728E+05, 0.12217E+05,
     + 0.12722E+05, 0.13242E+05, 0.13778E+05, 0.14331E+05, 0.14900E+05,
     + 0.15486E+05, 0.16091E+05, 0.16713E+05, 0.17353E+05, 0.18012E+05,
     + 0.18691E+05, 0.19389E+05, 0.20108E+05, 0.20847E+05, 0.21607E+05,
     + 0.22388E+05, 0.23191E+05, 0.24017E+05, 0.24866E+05, 0.25738E+05,
     + 0.26633E+05, 0.27553E+05, 0.28498E+05, 0.29468E+05, 0.30464E+05,
     + 0.31486E+05, 0.32536E+05, 0.33612E+05, 0.34716E+05, 0.35849E+05,
     + 0.37011E+05, 0.38202E+05, 0.39424E+05, 0.40676E+05, 0.41959E+05,
     + 0.43274E+05, 0.44622E+05, 0.46002E+05, 0.47416E+05, 0.48864E+05,
     + 0.50348E+05, 0.51866E+05, 0.53421E+05, 0.55012E+05, 0.56640E+05,
     + 0.58307E+05, 0.60012E+05, 0.61757E+05, 0.63541E+05, 0.65366E+05,
     + 0.67233E+05, 0.69141E+05, 0.71092E+05, 0.73087E+05, 0.75125E+05,
     + 0.77209E+05, 0.79338E+05, 0.81513E+05, 0.83736E+05, 0.86006E+05,
     + 0.88324E+05, 0.90693E+05, 0.93111E+05, 0.95580E+05, 0.98100E+05,
     + 0.10067E+06/
c...        --       172
      data (QofT( 6,J),J=1,119)/ 0.49379E+03, 0.82021E+03, 0.11980E+04,
     + 0.16211E+04, 0.20851E+04, 0.25870E+04, 0.31242E+04, 0.36949E+04,
     + 0.42977E+04, 0.49317E+04, 0.55963E+04, 0.62911E+04, 0.70164E+04,
     + 0.77722E+04, 0.85591E+04, 0.93777E+04, 0.10228E+05, 0.11112E+05,
     + 0.12030E+05, 0.12983E+05, 0.13971E+05, 0.14997E+05, 0.16061E+05,
     + 0.17163E+05, 0.18306E+05, 0.19491E+05, 0.20719E+05, 0.21991E+05,
     + 0.23309E+05, 0.24673E+05, 0.26086E+05, 0.27549E+05, 0.29064E+05,
     + 0.30631E+05, 0.32254E+05, 0.33932E+05, 0.35669E+05, 0.37464E+05,
     + 0.39321E+05, 0.41242E+05, 0.43227E+05, 0.45279E+05, 0.47399E+05,
     + 0.49589E+05, 0.51852E+05, 0.54189E+05, 0.56602E+05, 0.59094E+05,
     + 0.61666E+05, 0.64320E+05, 0.67058E+05, 0.69883E+05, 0.72796E+05,
     + 0.75801E+05, 0.78899E+05, 0.82092E+05, 0.85382E+05, 0.88773E+05,
     + 0.92266E+05, 0.95863E+05, 0.99568E+05, 0.10338E+06, 0.10731E+06,
     + 0.11135E+06, 0.11551E+06, 0.11979E+06, 0.12419E+06, 0.12871E+06,
     + 0.13337E+06, 0.13815E+06, 0.14307E+06, 0.14812E+06, 0.15331E+06,
     + 0.15865E+06, 0.16412E+06, 0.16975E+06, 0.17553E+06, 0.18146E+06,
     + 0.18754E+06, 0.19379E+06, 0.20020E+06, 0.20678E+06, 0.21352E+06,
     + 0.22044E+06, 0.22753E+06, 0.23480E+06, 0.24226E+06, 0.24990E+06,
     + 0.25773E+06, 0.26575E+06, 0.27397E+06, 0.28239E+06, 0.29102E+06,
     + 0.29985E+06, 0.30889E+06, 0.31814E+06, 0.32762E+06, 0.33731E+06,
     + 0.34724E+06, 0.35739E+06, 0.36777E+06, 0.37840E+06, 0.38926E+06,
     + 0.40038E+06, 0.41174E+06, 0.42335E+06, 0.43523E+06, 0.44737E+06,
     + 0.45977E+06, 0.47245E+06, 0.48540E+06, 0.49863E+06, 0.51214E+06,
     + 0.52595E+06, 0.54005E+06, 0.55444E+06, 0.56914E+06, 0.58415E+06,
     + 0.59947E+06/
  
      eps=0.01
c
      gsi = xgj(iso)
      do I=1,NT
         Q(I)=QofT(iso,I)
      end do
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
         Qt = -1.
         write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
         go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_CO2   (                       
     I T,       ! temperature in K 
     I iso,       ! isotope code (HITRAN INDEX)
     O gsi,       ! state independent nuclear degeneracyfactor
     O QT)       ! Total Internal Partition Function
 
      implicit DOUBLE PRECISION (a-h,o-z)
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj(11), QofT(11,119),Q(NT)
      data xgj/ 1.,2.,1.,6.,2.,12.,1.,6.,1.,2.,12./
c...      CO2
c...        --       626
      data (QofT( 1,J),J=1,119)/ 0.53642E+02, 0.75947E+02, 0.98292E+02,
     + 0.12078E+03, 0.14364E+03, 0.16714E+03, 0.19160E+03, 0.21731E+03,
     + 0.24454E+03, 0.27355E+03, 0.30456E+03, 0.33778E+03, 0.37343E+03,
     + 0.41170E+03, 0.45280E+03, 0.49692E+03, 0.54427E+03, 0.59505E+03,
     + 0.64948E+03, 0.70779E+03, 0.77019E+03, 0.83693E+03, 0.90825E+03,
     + 0.98440E+03, 0.10656E+04, 0.11522E+04, 0.12445E+04, 0.13427E+04,
     + 0.14471E+04, 0.15580E+04, 0.16759E+04, 0.18009E+04, 0.19334E+04,
     + 0.20739E+04, 0.22225E+04, 0.23798E+04, 0.25462E+04, 0.27219E+04,
     + 0.29074E+04, 0.31032E+04, 0.33097E+04, 0.35272E+04, 0.37564E+04,
     + 0.39976E+04, 0.42514E+04, 0.45181E+04, 0.47985E+04, 0.50929E+04,
     + 0.54019E+04, 0.57260E+04, 0.60659E+04, 0.64221E+04, 0.67952E+04,
     + 0.71859E+04, 0.75946E+04, 0.80222E+04, 0.84691E+04, 0.89362E+04,
     + 0.94241E+04, 0.99335E+04, 0.10465E+05, 0.11020E+05, 0.11598E+05,
     + 0.12201E+05, 0.12828E+05, 0.13482E+05, 0.14163E+05, 0.14872E+05,
     + 0.15609E+05, 0.16376E+05, 0.17173E+05, 0.18001E+05, 0.18861E+05,
     + 0.19754E+05, 0.20682E+05, 0.21644E+05, 0.22643E+05, 0.23678E+05,
     + 0.24752E+05, 0.25865E+05, 0.27018E+05, 0.28212E+05, 0.29449E+05,
     + 0.30730E+05, 0.32055E+05, 0.33426E+05, 0.34845E+05, 0.36312E+05,
     + 0.37828E+05, 0.39395E+05, 0.41015E+05, 0.42688E+05, 0.44416E+05,
     + 0.46199E+05, 0.48041E+05, 0.49942E+05, 0.51902E+05, 0.53925E+05,
     + 0.56011E+05, 0.58162E+05, 0.60379E+05, 0.62664E+05, 0.65019E+05,
     + 0.67444E+05, 0.69942E+05, 0.72515E+05, 0.75163E+05, 0.77890E+05,
     + 0.80695E+05, 0.83582E+05, 0.86551E+05, 0.89605E+05, 0.92746E+05,
     + 0.95975E+05, 0.99294E+05, 0.10271E+06, 0.10621E+06, 0.10981E+06,
     + 0.11351E+06/
c...        --       636
      data (QofT( 2,J),J=1,119)/ 0.10728E+03, 0.15189E+03, 0.19659E+03,
     + 0.24164E+03, 0.28753E+03, 0.33486E+03, 0.38429E+03, 0.43643E+03,
     + 0.49184E+03, 0.55104E+03, 0.61449E+03, 0.68263E+03, 0.75589E+03,
     + 0.83468E+03, 0.91943E+03, 0.10106E+04, 0.11085E+04, 0.12137E+04,
     + 0.13266E+04, 0.14477E+04, 0.15774E+04, 0.17163E+04, 0.18649E+04,
     + 0.20237E+04, 0.21933E+04, 0.23743E+04, 0.25673E+04, 0.27729E+04,
     + 0.29917E+04, 0.32245E+04, 0.34718E+04, 0.37345E+04, 0.40132E+04,
     + 0.43087E+04, 0.46218E+04, 0.49533E+04, 0.53041E+04, 0.56749E+04,
     + 0.60668E+04, 0.64805E+04, 0.69171E+04, 0.73774E+04, 0.78626E+04,
     + 0.83736E+04, 0.89114E+04, 0.94772E+04, 0.10072E+05, 0.10697E+05,
     + 0.11353E+05, 0.12042E+05, 0.12765E+05, 0.13523E+05, 0.14317E+05,
     + 0.15148E+05, 0.16019E+05, 0.16930E+05, 0.17883E+05, 0.18879E+05,
     + 0.19920E+05, 0.21008E+05, 0.22143E+05, 0.23328E+05, 0.24563E+05,
     + 0.25852E+05, 0.27195E+05, 0.28594E+05, 0.30051E+05, 0.31568E+05,
     + 0.33146E+05, 0.34788E+05, 0.36496E+05, 0.38271E+05, 0.40115E+05,
     + 0.42031E+05, 0.44021E+05, 0.46086E+05, 0.48230E+05, 0.50453E+05,
     + 0.52759E+05, 0.55150E+05, 0.57628E+05, 0.60195E+05, 0.62854E+05,
     + 0.65608E+05, 0.68459E+05, 0.71409E+05, 0.74461E+05, 0.77618E+05,
     + 0.80883E+05, 0.84258E+05, 0.87746E+05, 0.91350E+05, 0.95073E+05,
     + 0.98918E+05, 0.10289E+06, 0.10698E+06, 0.11121E+06, 0.11558E+06,
     + 0.12008E+06, 0.12472E+06, 0.12950E+06, 0.13443E+06, 0.13952E+06,
     + 0.14475E+06, 0.15015E+06, 0.15571E+06, 0.16143E+06, 0.16732E+06,
     + 0.17338E+06, 0.17962E+06, 0.18604E+06, 0.19264E+06, 0.19943E+06,
     + 0.20642E+06, 0.21360E+06, 0.22098E+06, 0.22856E+06, 0.23636E+06,
     + 0.24436E+06/
c...        --       628
      data (QofT( 3,J),J=1,119)/ 0.11368E+03, 0.16096E+03, 0.20833E+03,
     + 0.25603E+03, 0.30452E+03, 0.35442E+03, 0.40640E+03, 0.46110E+03,
     + 0.51910E+03, 0.58093E+03, 0.64709E+03, 0.71804E+03, 0.79422E+03,
     + 0.87607E+03, 0.96402E+03, 0.10585E+04, 0.11600E+04, 0.12689E+04,
     + 0.13857E+04, 0.15108E+04, 0.16449E+04, 0.17883E+04, 0.19416E+04,
     + 0.21054E+04, 0.22803E+04, 0.24668E+04, 0.26655E+04, 0.28770E+04,
     + 0.31021E+04, 0.33414E+04, 0.35956E+04, 0.38654E+04, 0.41516E+04,
     + 0.44549E+04, 0.47761E+04, 0.51160E+04, 0.54755E+04, 0.58555E+04,
     + 0.62568E+04, 0.66804E+04, 0.71273E+04, 0.75982E+04, 0.80944E+04,
     + 0.86169E+04, 0.91666E+04, 0.97446E+04, 0.10352E+05, 0.10990E+05,
     + 0.11660E+05, 0.12363E+05, 0.13101E+05, 0.13874E+05, 0.14683E+05,
     + 0.15531E+05, 0.16418E+05, 0.17347E+05, 0.18317E+05, 0.19332E+05,
     + 0.20392E+05, 0.21499E+05, 0.22654E+05, 0.23859E+05, 0.25116E+05,
     + 0.26426E+05, 0.27792E+05, 0.29214E+05, 0.30695E+05, 0.32236E+05,
     + 0.33840E+05, 0.35508E+05, 0.37242E+05, 0.39045E+05, 0.40917E+05,
     + 0.42862E+05, 0.44881E+05, 0.46977E+05, 0.49152E+05, 0.51407E+05,
     + 0.53746E+05, 0.56171E+05, 0.58683E+05, 0.61286E+05, 0.63981E+05,
     + 0.66772E+05, 0.69661E+05, 0.72650E+05, 0.75742E+05, 0.78940E+05,
     + 0.82246E+05, 0.85664E+05, 0.89196E+05, 0.92845E+05, 0.96613E+05,
     + 0.10050E+06, 0.10452E+06, 0.10867E+06, 0.11295E+06, 0.11736E+06,
     + 0.12191E+06, 0.12661E+06, 0.13145E+06, 0.13643E+06, 0.14157E+06,
     + 0.14687E+06, 0.15232E+06, 0.15794E+06, 0.16372E+06, 0.16968E+06,
     + 0.17580E+06, 0.18211E+06, 0.18859E+06, 0.19526E+06, 0.20213E+06,
     + 0.20918E+06, 0.21643E+06, 0.22388E+06, 0.23154E+06, 0.23941E+06,
     + 0.24750E+06/
c...        --       627
      data (QofT( 4,J),J=1,119)/ 0.66338E+03, 0.93923E+03, 0.12156E+04,
     + 0.14938E+04, 0.17766E+04, 0.20676E+04, 0.23705E+04, 0.26891E+04,
     + 0.30267E+04, 0.33866E+04, 0.37714E+04, 0.41839E+04, 0.46267E+04,
     + 0.51023E+04, 0.56132E+04, 0.61618E+04, 0.67508E+04, 0.73827E+04,
     + 0.80603E+04, 0.87863E+04, 0.95636E+04, 0.10395E+05, 0.11284E+05,
     + 0.12233E+05, 0.13246E+05, 0.14326E+05, 0.15477E+05, 0.16702E+05,
     + 0.18005E+05, 0.19390E+05, 0.20861E+05, 0.22422E+05, 0.24077E+05,
     + 0.25832E+05, 0.27689E+05, 0.29655E+05, 0.31734E+05, 0.33931E+05,
     + 0.36250E+05, 0.38698E+05, 0.41280E+05, 0.44002E+05, 0.46869E+05,
     + 0.49886E+05, 0.53062E+05, 0.56400E+05, 0.59909E+05, 0.63594E+05,
     + 0.67462E+05, 0.71521E+05, 0.75777E+05, 0.80238E+05, 0.84911E+05,
     + 0.89804E+05, 0.94925E+05, 0.10028E+06, 0.10588E+06, 0.11173E+06,
     + 0.11785E+06, 0.12423E+06, 0.13090E+06, 0.13785E+06, 0.14510E+06,
     + 0.15265E+06, 0.16053E+06, 0.16873E+06, 0.17727E+06, 0.18615E+06,
     + 0.19540E+06, 0.20501E+06, 0.21501E+06, 0.22540E+06, 0.23619E+06,
     + 0.24740E+06, 0.25904E+06, 0.27112E+06, 0.28365E+06, 0.29664E+06,
     + 0.31012E+06, 0.32409E+06, 0.33856E+06, 0.35356E+06, 0.36908E+06,
     + 0.38516E+06, 0.40180E+06, 0.41902E+06, 0.43683E+06, 0.45525E+06,
     + 0.47429E+06, 0.49397E+06, 0.51431E+06, 0.53532E+06, 0.55702E+06,
     + 0.57943E+06, 0.60256E+06, 0.62644E+06, 0.65107E+06, 0.67648E+06,
     + 0.70269E+06, 0.72972E+06, 0.75758E+06, 0.78629E+06, 0.81588E+06,
     + 0.84636E+06, 0.87775E+06, 0.91008E+06, 0.94337E+06, 0.97763E+06,
     + 0.10129E+07, 0.10492E+07, 0.10865E+07, 0.11249E+07, 0.11644E+07,
     + 0.12050E+07, 0.12467E+07, 0.12896E+07, 0.13337E+07, 0.13789E+07,
     + 0.14255E+07/
c...        --       638
      data (QofT( 5,J),J=1,119)/ 0.22737E+03, 0.32194E+03, 0.41671E+03,
     + 0.51226E+03, 0.60963E+03, 0.71017E+03, 0.81528E+03, 0.92628E+03,
     + 0.10444E+04, 0.11707E+04, 0.13061E+04, 0.14518E+04, 0.16085E+04,
     + 0.17772E+04, 0.19588E+04, 0.21542E+04, 0.23644E+04, 0.25903E+04,
     + 0.28330E+04, 0.30934E+04, 0.33726E+04, 0.36717E+04, 0.39918E+04,
     + 0.43342E+04, 0.47001E+04, 0.50907E+04, 0.55074E+04, 0.59515E+04,
     + 0.64244E+04, 0.69276E+04, 0.74626E+04, 0.80310E+04, 0.86344E+04,
     + 0.92744E+04, 0.99528E+04, 0.10671E+05, 0.11432E+05, 0.12236E+05,
     + 0.13086E+05, 0.13984E+05, 0.14932E+05, 0.15932E+05, 0.16985E+05,
     + 0.18096E+05, 0.19265E+05, 0.20495E+05, 0.21788E+05, 0.23148E+05,
     + 0.24576E+05, 0.26075E+05, 0.27648E+05, 0.29298E+05, 0.31027E+05,
     + 0.32839E+05, 0.34736E+05, 0.36721E+05, 0.38798E+05, 0.40970E+05,
     + 0.43240E+05, 0.45611E+05, 0.48087E+05, 0.50671E+05, 0.53368E+05,
     + 0.56180E+05, 0.59111E+05, 0.62165E+05, 0.65347E+05, 0.68659E+05,
     + 0.72107E+05, 0.75694E+05, 0.79425E+05, 0.83303E+05, 0.87334E+05,
     + 0.91522E+05, 0.95872E+05, 0.10039E+06, 0.10507E+06, 0.10994E+06,
     + 0.11498E+06, 0.12021E+06, 0.12563E+06, 0.13125E+06, 0.13707E+06,
     + 0.14309E+06, 0.14933E+06, 0.15579E+06, 0.16247E+06, 0.16938E+06,
     + 0.17653E+06, 0.18392E+06, 0.19156E+06, 0.19946E+06, 0.20761E+06,
     + 0.21604E+06, 0.22473E+06, 0.23371E+06, 0.24298E+06, 0.25254E+06,
     + 0.26240E+06, 0.27258E+06, 0.28307E+06, 0.29388E+06, 0.30502E+06,
     + 0.31651E+06, 0.32834E+06, 0.34052E+06, 0.35307E+06, 0.36599E+06,
     + 0.37929E+06, 0.39298E+06, 0.40706E+06, 0.42155E+06, 0.43645E+06,
     + 0.45178E+06, 0.46753E+06, 0.48373E+06, 0.50038E+06, 0.51748E+06,
     + 0.53506E+06/
c...        --       637
      data (QofT( 6,J),J=1,119)/ 0.13267E+04, 0.18785E+04, 0.24314E+04,
     + 0.29888E+04, 0.35566E+04, 0.41426E+04, 0.47550E+04, 0.54013E+04,
     + 0.60886E+04, 0.68232E+04, 0.76109E+04, 0.84574E+04, 0.93678E+04,
     + 0.10348E+05, 0.11402E+05, 0.12536E+05, 0.13755E+05, 0.15065E+05,
     + 0.16471E+05, 0.17980E+05, 0.19598E+05, 0.21330E+05, 0.23184E+05,
     + 0.25166E+05, 0.27283E+05, 0.29543E+05, 0.31953E+05, 0.34521E+05,
     + 0.37256E+05, 0.40164E+05, 0.43256E+05, 0.46541E+05, 0.50026E+05,
     + 0.53723E+05, 0.57641E+05, 0.61790E+05, 0.66180E+05, 0.70823E+05,
     + 0.75729E+05, 0.80910E+05, 0.86378E+05, 0.92145E+05, 0.98224E+05,
     + 0.10463E+06, 0.11137E+06, 0.11846E+06, 0.12592E+06, 0.13375E+06,
     + 0.14198E+06, 0.15062E+06, 0.15969E+06, 0.16920E+06, 0.17916E+06,
     + 0.18959E+06, 0.20052E+06, 0.21196E+06, 0.22392E+06, 0.23642E+06,
     + 0.24949E+06, 0.26314E+06, 0.27740E+06, 0.29227E+06, 0.30779E+06,
     + 0.32398E+06, 0.34085E+06, 0.35842E+06, 0.37673E+06, 0.39579E+06,
     + 0.41563E+06, 0.43626E+06, 0.45772E+06, 0.48003E+06, 0.50322E+06,
     + 0.52730E+06, 0.55232E+06, 0.57829E+06, 0.60524E+06, 0.63320E+06,
     + 0.66219E+06, 0.69226E+06, 0.72342E+06, 0.75571E+06, 0.78916E+06,
     + 0.82380E+06, 0.85966E+06, 0.89678E+06, 0.93518E+06, 0.97490E+06,
     + 0.10160E+07, 0.10585E+07, 0.11023E+07, 0.11477E+07, 0.11946E+07,
     + 0.12430E+07, 0.12929E+07, 0.13445E+07, 0.13977E+07, 0.14526E+07,
     + 0.15093E+07, 0.15677E+07, 0.16280E+07, 0.16901E+07, 0.17541E+07,
     + 0.18200E+07, 0.18880E+07, 0.19579E+07, 0.20300E+07, 0.21042E+07,
     + 0.21805E+07, 0.22591E+07, 0.23400E+07, 0.24232E+07, 0.25087E+07,
     + 0.25967E+07, 0.26871E+07, 0.27801E+07, 0.28757E+07, 0.29739E+07,
     + 0.30747E+07/
c...        --       828
      data (QofT( 7,J),J=1,119)/ 0.60334E+02, 0.85430E+02, 0.11058E+03,
     + 0.13590E+03, 0.16167E+03, 0.18821E+03, 0.21588E+03, 0.24502E+03,
     + 0.27595E+03, 0.30896E+03, 0.34431E+03, 0.38225E+03, 0.42301E+03,
     + 0.46684E+03, 0.51397E+03, 0.56464E+03, 0.61907E+03, 0.67753E+03,
     + 0.74027E+03, 0.80753E+03, 0.87961E+03, 0.95676E+03, 0.10393E+04,
     + 0.11275E+04, 0.12217E+04, 0.13222E+04, 0.14293E+04, 0.15434E+04,
     + 0.16648E+04, 0.17940E+04, 0.19312E+04, 0.20769E+04, 0.22315E+04,
     + 0.23954E+04, 0.25691E+04, 0.27529E+04, 0.29474E+04, 0.31530E+04,
     + 0.33702E+04, 0.35995E+04, 0.38414E+04, 0.40965E+04, 0.43654E+04,
     + 0.46484E+04, 0.49464E+04, 0.52598E+04, 0.55892E+04, 0.59353E+04,
     + 0.62988E+04, 0.66803E+04, 0.70804E+04, 0.74998E+04, 0.79394E+04,
     + 0.83998E+04, 0.88817E+04, 0.93859E+04, 0.99132E+04, 0.10464E+05,
     + 0.11040E+05, 0.11642E+05, 0.12270E+05, 0.12925E+05, 0.13609E+05,
     + 0.14321E+05, 0.15064E+05, 0.15838E+05, 0.16643E+05, 0.17482E+05,
     + 0.18355E+05, 0.19263E+05, 0.20207E+05, 0.21188E+05, 0.22208E+05,
     + 0.23267E+05, 0.24366E+05, 0.25508E+05, 0.26692E+05, 0.27921E+05,
     + 0.29195E+05, 0.30516E+05, 0.31886E+05, 0.33304E+05, 0.34773E+05,
     + 0.36294E+05, 0.37869E+05, 0.39499E+05, 0.41185E+05, 0.42929E+05,
     + 0.44732E+05, 0.46596E+05, 0.48522E+05, 0.50513E+05, 0.52569E+05,
     + 0.54692E+05, 0.56884E+05, 0.59146E+05, 0.61481E+05, 0.63890E+05,
     + 0.66375E+05, 0.68937E+05, 0.71578E+05, 0.74301E+05, 0.77107E+05,
     + 0.79998E+05, 0.82976E+05, 0.86043E+05, 0.89201E+05, 0.92452E+05,
     + 0.95799E+05, 0.99242E+05, 0.10278E+06, 0.10643E+06, 0.11018E+06,
     + 0.11403E+06, 0.11799E+06, 0.12206E+06, 0.12625E+06, 0.13055E+06,
     + 0.13497E+06/
c...        --       728
      data (QofT( 8,J),J=1,119)/ 0.70354E+03, 0.99615E+03, 0.12893E+04,
     + 0.15846E+04, 0.18848E+04, 0.21940E+04, 0.25162E+04, 0.28554E+04,
     + 0.32152E+04, 0.35991E+04, 0.40099E+04, 0.44507E+04, 0.49242E+04,
     + 0.54332E+04, 0.59802E+04, 0.65681E+04, 0.71996E+04, 0.78776E+04,
     + 0.86050E+04, 0.93847E+04, 0.10220E+05, 0.11114E+05, 0.12070E+05,
     + 0.13091E+05, 0.14182E+05, 0.15345E+05, 0.16585E+05, 0.17906E+05,
     + 0.19311E+05, 0.20805E+05, 0.22393E+05, 0.24078E+05, 0.25865E+05,
     + 0.27760E+05, 0.29768E+05, 0.31893E+05, 0.34140E+05, 0.36516E+05,
     + 0.39025E+05, 0.41674E+05, 0.44469E+05, 0.47416E+05, 0.50520E+05,
     + 0.53789E+05, 0.57229E+05, 0.60847E+05, 0.64650E+05, 0.68645E+05,
     + 0.72840E+05, 0.77242E+05, 0.81859E+05, 0.86699E+05, 0.91770E+05,
     + 0.97081E+05, 0.10264E+06, 0.10846E+06, 0.11454E+06, 0.12090E+06,
     + 0.12754E+06, 0.13447E+06, 0.14171E+06, 0.14927E+06, 0.15715E+06,
     + 0.16536E+06, 0.17392E+06, 0.18284E+06, 0.19213E+06, 0.20179E+06,
     + 0.21185E+06, 0.22231E+06, 0.23319E+06, 0.24450E+06, 0.25625E+06,
     + 0.26845E+06, 0.28112E+06, 0.29427E+06, 0.30791E+06, 0.32206E+06,
     + 0.33674E+06, 0.35196E+06, 0.36772E+06, 0.38406E+06, 0.40098E+06,
     + 0.41850E+06, 0.43663E+06, 0.45539E+06, 0.47480E+06, 0.49488E+06,
     + 0.51564E+06, 0.53710E+06, 0.55928E+06, 0.58219E+06, 0.60586E+06,
     + 0.63029E+06, 0.65553E+06, 0.68157E+06, 0.70844E+06, 0.73616E+06,
     + 0.76476E+06, 0.79424E+06, 0.82464E+06, 0.85597E+06, 0.88826E+06,
     + 0.92153E+06, 0.95580E+06, 0.99108E+06, 0.10274E+07, 0.10648E+07,
     + 0.11033E+07, 0.11429E+07, 0.11837E+07, 0.12256E+07, 0.12687E+07,
     + 0.13131E+07, 0.13586E+07, 0.14055E+07, 0.14536E+07, 0.15031E+07,
     + 0.15539E+07/
c...        --       727
      data (QofT( 9,J),J=1,119)/ 0.20518E+04, 0.29051E+04, 0.37601E+04,
     + 0.46209E+04, 0.54961E+04, 0.63969E+04, 0.73353E+04, 0.83227E+04,
     + 0.93698E+04, 0.10486E+05, 0.11681E+05, 0.12962E+05, 0.14337E+05,
     + 0.15815E+05, 0.17403E+05, 0.19110E+05, 0.20942E+05, 0.22909E+05,
     + 0.25018E+05, 0.27278E+05, 0.29699E+05, 0.32290E+05, 0.35060E+05,
     + 0.38019E+05, 0.41177E+05, 0.44545E+05, 0.48135E+05, 0.51957E+05,
     + 0.56023E+05, 0.60346E+05, 0.64938E+05, 0.69812E+05, 0.74981E+05,
     + 0.80461E+05, 0.86264E+05, 0.92406E+05, 0.98902E+05, 0.10577E+06,
     + 0.11302E+06, 0.12067E+06, 0.12875E+06, 0.13726E+06, 0.14622E+06,
     + 0.15566E+06, 0.16559E+06, 0.17604E+06, 0.18702E+06, 0.19855E+06,
     + 0.21066E+06, 0.22336E+06, 0.23669E+06, 0.25065E+06, 0.26528E+06,
     + 0.28061E+06, 0.29664E+06, 0.31342E+06, 0.33096E+06, 0.34930E+06,
     + 0.36845E+06, 0.38845E+06, 0.40933E+06, 0.43111E+06, 0.45383E+06,
     + 0.47751E+06, 0.50219E+06, 0.52790E+06, 0.55466E+06, 0.58252E+06,
     + 0.61151E+06, 0.64166E+06, 0.67300E+06, 0.70558E+06, 0.73943E+06,
     + 0.77458E+06, 0.81108E+06, 0.84896E+06, 0.88827E+06, 0.92904E+06,
     + 0.97131E+06, 0.10151E+07, 0.10605E+07, 0.11076E+07, 0.11563E+07,
     + 0.12068E+07, 0.12590E+07, 0.13130E+07, 0.13689E+07, 0.14267E+07,
     + 0.14865E+07, 0.15483E+07, 0.16121E+07, 0.16781E+07, 0.17462E+07,
     + 0.18165E+07, 0.18892E+07, 0.19641E+07, 0.20415E+07, 0.21213E+07,
     + 0.22036E+07, 0.22884E+07, 0.23759E+07, 0.24661E+07, 0.25590E+07,
     + 0.26547E+07, 0.27533E+07, 0.28549E+07, 0.29594E+07, 0.30670E+07,
     + 0.31778E+07, 0.32918E+07, 0.34090E+07, 0.35296E+07, 0.36536E+07,
     + 0.37812E+07, 0.39123E+07, 0.40470E+07, 0.41855E+07, 0.43278E+07,
     + 0.44739E+07/
c...        --       838
      data (QofT(10,J),J=1,119)/ 0.12066E+03, 0.17085E+03, 0.22116E+03,
     + 0.27190E+03, 0.32364E+03, 0.37711E+03, 0.43305E+03, 0.49219E+03,
     + 0.55516E+03, 0.62256E+03, 0.69492E+03, 0.77276E+03, 0.85657E+03,
     + 0.94685E+03, 0.10441E+04, 0.11488E+04, 0.12614E+04, 0.13826E+04,
     + 0.15127E+04, 0.16525E+04, 0.18024E+04, 0.19630E+04, 0.21351E+04,
     + 0.23191E+04, 0.25158E+04, 0.27260E+04, 0.29502E+04, 0.31892E+04,
     + 0.34438E+04, 0.37148E+04, 0.40031E+04, 0.43094E+04, 0.46346E+04,
     + 0.49797E+04, 0.53455E+04, 0.57331E+04, 0.61434E+04, 0.65775E+04,
     + 0.70364E+04, 0.75212E+04, 0.80330E+04, 0.85730E+04, 0.91424E+04,
     + 0.97423E+04, 0.10374E+05, 0.11039E+05, 0.11738E+05, 0.12474E+05,
     + 0.13246E+05, 0.14057E+05, 0.14908E+05, 0.15801E+05, 0.16737E+05,
     + 0.17717E+05, 0.18744E+05, 0.19819E+05, 0.20944E+05, 0.22120E+05,
     + 0.23349E+05, 0.24634E+05, 0.25975E+05, 0.27376E+05, 0.28837E+05,
     + 0.30361E+05, 0.31950E+05, 0.33605E+05, 0.35330E+05, 0.37126E+05,
     + 0.38996E+05, 0.40942E+05, 0.42965E+05, 0.45069E+05, 0.47256E+05,
     + 0.49528E+05, 0.51888E+05, 0.54338E+05, 0.56882E+05, 0.59521E+05,
     + 0.62259E+05, 0.65097E+05, 0.68040E+05, 0.71090E+05, 0.74249E+05,
     + 0.77522E+05, 0.80910E+05, 0.84417E+05, 0.88046E+05, 0.91801E+05,
     + 0.95684E+05, 0.99699E+05, 0.10385E+06, 0.10814E+06, 0.11257E+06,
     + 0.11715E+06, 0.12187E+06, 0.12675E+06, 0.13179E+06, 0.13699E+06,
     + 0.14235E+06, 0.14788E+06, 0.15358E+06, 0.15946E+06, 0.16552E+06,
     + 0.17176E+06, 0.17819E+06, 0.18482E+06, 0.19164E+06, 0.19867E+06,
     + 0.20590E+06, 0.21335E+06, 0.22101E+06, 0.22889E+06, 0.23699E+06,
     + 0.24533E+06, 0.25390E+06, 0.26271E+06, 0.27177E+06, 0.28108E+06,
     + 0.29064E+06/
c...        --       837
      data (QofT(11,J),J=1,119)/ 0.14071E+04, 0.19923E+04, 0.25789E+04,
     + 0.31704E+04, 0.37733E+04, 0.43962E+04, 0.50477E+04, 0.57360E+04,
     + 0.64687E+04, 0.72525E+04, 0.80938E+04, 0.89984E+04, 0.99723E+04,
     + 0.11021E+05, 0.12150E+05, 0.13366E+05, 0.14673E+05, 0.16079E+05,
     + 0.17589E+05, 0.19211E+05, 0.20949E+05, 0.22812E+05, 0.24807E+05,
     + 0.26940E+05, 0.29221E+05, 0.31656E+05, 0.34254E+05, 0.37023E+05,
     + 0.39972E+05, 0.43111E+05, 0.46449E+05, 0.49996E+05, 0.53762E+05,
     + 0.57756E+05, 0.61991E+05, 0.66477E+05, 0.71226E+05, 0.76249E+05,
     + 0.81558E+05, 0.87167E+05, 0.93088E+05, 0.99334E+05, 0.10592E+06,
     + 0.11286E+06, 0.12016E+06, 0.12785E+06, 0.13594E+06, 0.14444E+06,
     + 0.15337E+06, 0.16274E+06, 0.17258E+06, 0.18290E+06, 0.19371E+06,
     + 0.20504E+06, 0.21691E+06, 0.22933E+06, 0.24233E+06, 0.25592E+06,
     + 0.27012E+06, 0.28496E+06, 0.30046E+06, 0.31663E+06, 0.33351E+06,
     + 0.35111E+06, 0.36946E+06, 0.38858E+06, 0.40850E+06, 0.42924E+06,
     + 0.45083E+06, 0.47329E+06, 0.49666E+06, 0.52095E+06, 0.54620E+06,
     + 0.57243E+06, 0.59967E+06, 0.62796E+06, 0.65732E+06, 0.68778E+06,
     + 0.71938E+06, 0.75214E+06, 0.78611E+06, 0.82131E+06, 0.85777E+06,
     + 0.89553E+06, 0.93463E+06, 0.97511E+06, 0.10170E+07, 0.10603E+07,
     + 0.11051E+07, 0.11514E+07, 0.11993E+07, 0.12488E+07, 0.12999E+07,
     + 0.13527E+07, 0.14073E+07, 0.14636E+07, 0.15217E+07, 0.15816E+07,
     + 0.16435E+07, 0.17072E+07, 0.17730E+07, 0.18408E+07, 0.19107E+07,
     + 0.19827E+07, 0.20569E+07, 0.21334E+07, 0.22121E+07, 0.22931E+07,
     + 0.23765E+07, 0.24624E+07, 0.25507E+07, 0.26416E+07, 0.27351E+07,
     + 0.28312E+07, 0.29301E+07, 0.30317E+07, 0.31361E+07, 0.32434E+07,
     + 0.33537E+07/
  
      eps=0.01
c
      gsi = xgj(iso)
      do I=1,NT
         Q(I)=QofT(iso,I)
      end do
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
         Qt = -1.
         write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
         go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c     *****************
      Subroutine QT_CH4   (                       
     I T,       ! temperature in K 
     I iso,       ! isotope code (HITRAN INDEX)
     O gsi,       ! state independent nuclear degeneracyfactor
     O QT)       ! Total Internal Partition Function
 
      implicit DOUBLE PRECISION (a-h,o-z)
      parameter (NT=119)
      COMMON/Temperatures/tdat(NT)
      
      dimension xgj( 4), QofT( 4,119),Q(NT)
      data xgj/ 1.,2.,3.,6./
c...      CH4
c...        --       211
      data (QofT( 1,J),J=1,119)/ 0.54800E+02, 0.91500E+02, 0.13410E+03,
     + 0.18180E+03, 0.23410E+03, 0.29070E+03, 0.35140E+03, 0.41600E+03,
     + 0.48450E+03, 0.55720E+03, 0.63420E+03, 0.71600E+03, 0.80310E+03,
     + 0.89590E+03, 0.99520E+03, 0.11017E+04, 0.12161E+04, 0.13393E+04,
     + 0.14721E+04, 0.16155E+04, 0.17706E+04, 0.19384E+04, 0.21202E+04,
     + 0.23172E+04, 0.25307E+04, 0.27624E+04, 0.30137E+04, 0.32864E+04,
     + 0.35823E+04, 0.39034E+04, 0.42519E+04, 0.46300E+04, 0.50402E+04,
     + 0.54853E+04, 0.59679E+04, 0.64913E+04, 0.70588E+04, 0.76739E+04,
     + 0.83404E+04, 0.90625E+04, 0.98446E+04, 0.10691E+05, 0.11608E+05,
     + 0.12600E+05, 0.13674E+05, 0.14835E+05, 0.16090E+05, 0.17447E+05,
     + 0.18914E+05, 0.20500E+05, 0.22212E+05, 0.24063E+05, 0.26061E+05,
     + 0.28218E+05, 0.30548E+05, 0.33063E+05, 0.35778E+05, 0.38708E+05,
     + 0.41871E+05, 0.45284E+05, 0.48970E+05, 0.52940E+05, 0.57230E+05,
     + 0.61860E+05, 0.66860E+05, 0.72250E+05, 0.78070E+05, 0.84350E+05,
     + 0.91130E+05, 0.98450E+05, 0.10635E+06, 0.11488E+06, 0.12408E+06,
     + 0.13403E+06, 0.14480E+06, 0.15640E+06, 0.16890E+06, 0.18240E+06,
     + 0.19700E+06, 0.21280E+06, 0.22980E+06, 0.24830E+06, 0.26820E+06,
     + 0.28970E+06, 0.31290E+06, 0.33800E+06, 0.36520E+06, 0.39450E+06,
     + 0.42600E+06, 0.46000E+06, 0.49700E+06, 0.53700E+06, 0.58100E+06,
     + 0.62700E+06, 0.67800E+06, 0.73300E+06, 0.79200E+06, 0.85600E+06,
     + 0.92500E+06, 0.10000E+07, 0.10800E+07, 0.11670E+07, 0.12610E+07,
     + 0.13620E+07, 0.14720E+07, 0.15910E+07, 0.17190E+07, 0.18600E+07,
     + 0.20100E+07, 0.21700E+07, 0.23400E+07, 0.25300E+07, 0.27300E+07,
     + 0.29500E+07, 0.31800E+07, 0.34300E+07, 0.37000E+07, 0.39900E+07,
     + 0.42856E+07/
c...        --       311
      data (QofT( 2,J),J=1,119)/ 0.10958E+03, 0.18304E+03, 0.26818E+03,
     + 0.36356E+03, 0.46820E+03, 0.58141E+03, 0.70270E+03, 0.83186E+03,
     + 0.96893E+03, 0.11142E+04, 0.12682E+04, 0.14316E+04, 0.16055E+04,
     + 0.17909E+04, 0.19891E+04, 0.22016E+04, 0.24297E+04, 0.26752E+04,
     + 0.29399E+04, 0.32255E+04, 0.35342E+04, 0.38680E+04, 0.42294E+04,
     + 0.46208E+04, 0.50449E+04, 0.55046E+04, 0.60030E+04, 0.65434E+04,
     + 0.71293E+04, 0.77646E+04, 0.84535E+04, 0.92004E+04, 0.10010E+05,
     + 0.10888E+05, 0.11838E+05, 0.12869E+05, 0.13984E+05, 0.15193E+05,
     + 0.16501E+05, 0.17916E+05, 0.19448E+05, 0.21104E+05, 0.22895E+05,
     + 0.24830E+05, 0.26921E+05, 0.29180E+05, 0.31618E+05, 0.34250E+05,
     + 0.37090E+05, 0.40152E+05, 0.43454E+05, 0.47012E+05, 0.50845E+05,
     + 0.54973E+05, 0.59416E+05, 0.64197E+05, 0.69340E+05, 0.74870E+05,
     + 0.80813E+05, 0.87198E+05, 0.94055E+05, 0.10142E+06, 0.10932E+06,
     + 0.11779E+06, 0.12688E+06, 0.13662E+06, 0.14706E+06, 0.15824E+06,
     + 0.17021E+06, 0.18302E+06, 0.19673E+06, 0.21139E+06, 0.22706E+06,
     + 0.24381E+06, 0.26171E+06, 0.28082E+06, 0.30122E+06, 0.32299E+06,
     + 0.34621E+06, 0.37097E+06, 0.39737E+06, 0.42551E+06, 0.45548E+06,
     + 0.48739E+06, 0.52136E+06, 0.55752E+06, 0.59598E+06, 0.63688E+06,
     + 0.68036E+06, 0.72657E+06, 0.77566E+06, 0.82780E+06, 0.88316E+06,
     + 0.94191E+06, 0.10043E+07, 0.10704E+07, 0.11405E+07, 0.12148E+07,
     + 0.12936E+07, 0.13770E+07, 0.14654E+07, 0.15589E+07, 0.16579E+07,
     + 0.17627E+07, 0.18736E+07, 0.19908E+07, 0.21147E+07, 0.22456E+07,
     + 0.23840E+07, 0.25301E+07, 0.26844E+07, 0.28474E+07, 0.30193E+07,
     + 0.32007E+07, 0.33921E+07, 0.35939E+07, 0.38067E+07, 0.40310E+07,
     + 0.42673E+07/
c...        --       212
      data (QofT( 3,J),J=1,119)/ 0.44079E+03, 0.73786E+03, 0.10822E+04,
     + 0.14679E+04, 0.18913E+04, 0.23497E+04, 0.28415E+04, 0.33665E+04,
     + 0.39257E+04, 0.45211E+04, 0.51562E+04, 0.58349E+04, 0.65624E+04,
     + 0.73445E+04, 0.81872E+04, 0.90978E+04, 0.10084E+05, 0.11153E+05,
     + 0.12315E+05, 0.13579E+05, 0.14955E+05, 0.16455E+05, 0.18089E+05,
     + 0.19871E+05, 0.21816E+05, 0.23937E+05, 0.26251E+05, 0.28776E+05,
     + 0.31531E+05, 0.34535E+05, 0.37811E+05, 0.41384E+05, 0.45278E+05,
     + 0.49521E+05, 0.54144E+05, 0.59178E+05, 0.64657E+05, 0.70621E+05,
     + 0.77108E+05, 0.84161E+05, 0.91828E+05, 0.10016E+06, 0.10921E+06,
     + 0.11903E+06, 0.12968E+06, 0.14124E+06, 0.15378E+06, 0.16736E+06,
     + 0.18207E+06, 0.19800E+06, 0.21524E+06, 0.23389E+06, 0.25405E+06,
     + 0.27585E+06, 0.29939E+06, 0.32482E+06, 0.35226E+06, 0.38186E+06,
     + 0.41379E+06, 0.44821E+06, 0.48529E+06, 0.52522E+06, 0.56821E+06,
     + 0.61447E+06, 0.66422E+06, 0.71771E+06, 0.77519E+06, 0.83693E+06,
     + 0.90323E+06, 0.97438E+06, 0.10507E+07, 0.11326E+07, 0.12203E+07,
     + 0.13143E+07, 0.14150E+07, 0.15228E+07, 0.16382E+07, 0.17616E+07,
     + 0.18935E+07, 0.20346E+07, 0.21853E+07, 0.23463E+07, 0.25181E+07,
     + 0.27016E+07, 0.28973E+07, 0.31060E+07, 0.33284E+07, 0.35655E+07,
     + 0.38181E+07, 0.40870E+07, 0.43733E+07, 0.46780E+07, 0.50020E+07,
     + 0.53467E+07, 0.57130E+07, 0.61023E+07, 0.65158E+07, 0.69549E+07,
     + 0.74211E+07, 0.79158E+07, 0.84407E+07, 0.89973E+07, 0.95874E+07,
     + 0.10213E+08, 0.10875E+08, 0.11577E+08, 0.12320E+08, 0.13107E+08,
     + 0.13940E+08, 0.14820E+08, 0.15752E+08, 0.16736E+08, 0.17777E+08,
     + 0.18877E+08, 0.20038E+08, 0.21265E+08, 0.22560E+08, 0.23927E+08,
     + 0.25369E+08/
c...        --       312
      data (QofT( 4,J),J=1,119)/ 0.88231E+03, 0.14770E+04, 0.21661E+04,
     + 0.29384E+04, 0.37859E+04, 0.47034E+04, 0.56879E+04, 0.67388E+04,
     + 0.78581E+04, 0.90501E+04, 0.10321E+05, 0.11680E+05, 0.13136E+05,
     + 0.14702E+05, 0.16389E+05, 0.18212E+05, 0.20186E+05, 0.22328E+05,
     + 0.24654E+05, 0.27185E+05, 0.29941E+05, 0.32943E+05, 0.36216E+05,
     + 0.39786E+05, 0.43681E+05, 0.47930E+05, 0.52567E+05, 0.57625E+05,
     + 0.63144E+05, 0.69164E+05, 0.75730E+05, 0.82890E+05, 0.90693E+05,
     + 0.99198E+05, 0.10846E+06, 0.11855E+06, 0.12954E+06, 0.14149E+06,
     + 0.15450E+06, 0.16864E+06, 0.18402E+06, 0.20072E+06, 0.21886E+06,
     + 0.23856E+06, 0.25993E+06, 0.28312E+06, 0.30825E+06, 0.33550E+06,
     + 0.36501E+06, 0.39696E+06, 0.43155E+06, 0.46896E+06, 0.50942E+06,
     + 0.55315E+06, 0.60039E+06, 0.65141E+06, 0.70648E+06, 0.76589E+06,
     + 0.82997E+06, 0.89904E+06, 0.97346E+06, 0.10536E+07, 0.11399E+07,
     + 0.12327E+07, 0.13326E+07, 0.14400E+07, 0.15554E+07, 0.16793E+07,
     + 0.18124E+07, 0.19553E+07, 0.21085E+07, 0.22729E+07, 0.24490E+07,
     + 0.26378E+07, 0.28400E+07, 0.30565E+07, 0.32881E+07, 0.35360E+07,
     + 0.38010E+07, 0.40843E+07, 0.43870E+07, 0.47103E+07, 0.50555E+07,
     + 0.54239E+07, 0.58169E+07, 0.62361E+07, 0.66830E+07, 0.71592E+07,
     + 0.76666E+07, 0.82069E+07, 0.87820E+07, 0.93940E+07, 0.10045E+08,
     + 0.10737E+08, 0.11473E+08, 0.12256E+08, 0.13086E+08, 0.13969E+08,
     + 0.14905E+08, 0.15899E+08, 0.16954E+08, 0.18072E+08, 0.19258E+08,
     + 0.20515E+08, 0.21847E+08, 0.23257E+08, 0.24750E+08, 0.26331E+08,
     + 0.28004E+08, 0.29774E+08, 0.31646E+08, 0.33625E+08, 0.35716E+08,
     + 0.37926E+08, 0.40261E+08, 0.42726E+08, 0.45329E+08, 0.48077E+08,
     + 0.50975E+08/
  
      eps=0.01
c
      gsi = xgj(iso)
      do I=1,NT
         Q(I)=QofT(iso,I)
      end do
 
c
c...value depends on temperature range
      if(T.lt.70. .OR. T.gt.3000.) then
         Qt = -1.
         write(*,'(a)') '  OUT OF TEMPERATURE RANGE'
         go to 99
      endif
  
      call AtoB(T,Qt,Tdat,Q,NT)
c      
   99 return
      end
c
c***************************
      SUBROUTINE AtoB(aa,bb,A,B,npt)
c***************************
c...LaGrange 3- and 4-point interpolation
c...arrays A and B are the npt data points,  given aa, a value of the 
c...A variable, the routine will find the corresponding bb value
c
c...input:  aa
c...output: bb 
      implicit DOUBLE PRECISION (a-h,o-z)
      Parameter (Nmax=119)
      dimension A(Nmax),B(Nmax)
c
C 
c
      DO I=2,npt
         IF(A(I).GE.aa)THEN 
            IF(I.LT.3 .OR. I.EQ.npt) THEN
C     LaGrange three point interpolation 
               J = I
               IF(I.LT.3) J = 3
               IF(I.EQ.npT) J = npt
c.....do not divide by zero
               A0D1=A(J-2)-A(J-1)
               IF(abs(A0D1).LE.0.) A0D1=0.0001
               A0D2=A(J-2)-A(J)
               IF(abs(A0D2).LE.0.) A0D2=0.0001
               A1D1=A(J-1)-A(J-2)
               IF(abs(A1D1).LE.0.) A1D1=0.0001
               A1D2=A(J-1)-A(J)
               IF(abs(A1D2).LE.0.) A1D2=0.0001
               A2D1=A(J)-A(J-2)
               IF(abs(A2D1).LE.0.) A2D1=0.0001
               A2D2=A(J)-A(J-1)
               IF(abs(A2D2).LE.0.) A2D2=0.0001
c
               A0=(aa-A(J-1))*(aa-A(J))/(A0D1*A0D2)
               A1=(aa-A(J-2))*(aa-A(J))/(A1D1*A1D2)
               A2=(aa-A(J-2))*(aa-A(J-1))/(A2D1*A2D2)
c
               bb = A0*B(J-2) + A1*B(J-1) + A2*B(J)
c
            ELSE
C     LaGrange four point interpolation 
               J = I
c.....do not devide by zero
               A0D1=A(J-2)-A(J-1)
               IF(abs(A0D1).LE.0.) A0D1=0.0001
               A0D2=A(J-2)-A(J)
               IF(abs(A0D2).LE.0.) A0D2=0.0001
               A0D3 = (A(J-2)-A(J+1))
               IF(abs(A0D3).LE.0.) A0D3=0.0001
c
               A1D1=A(J-1)-A(J-2)
               IF(abs(A1D1).LE.0.) A1D1=0.0001
               A1D2=A(J-1)-A(J)
               IF(abs(A1D2).LE.0.) A1D2=0.0001
               A1D3 = A(J-1)-A(J+1)
               IF(abs(A1D3).LE.0.) A1D3=0.0001
c
               A2D1=A(J)-A(J-2)
               IF(abs(A2D1).LE.0.) A2D1=0.0001
               A2D2=A(J)-A(J-1)
               IF(abs(A2D2).LE.0.) A2D2=0.0001
               A2D3 = A(J)-A(J+1)
               IF(abs(A2D3).LE.0.) A2D3=0.0001
c
               A3D1 = A(J+1)-A(J-2)
               IF(abs(A3D1).LE.0.) A3D1=0.0001
               A3D2 = A(J+1)-A(J-1)
               IF(abs(A3D2).LE.0.) A3D2=0.0001
               A3D3 = A(J+1)-A(J)
               IF(abs(A3D3).LE.0.) A3D3=0.0001
c
               A0=(aa-A(J-1))*(aa-A(J))*(aa-A(J+1))
               A0=A0/(A0D1*A0D2*A0D3)
               A1=(aa-A(J-2))*(aa-A(J))*(aa-A(J+1))
               A1=A1/(A1D1*A1D2*A1D3)
               A2=(aa-A(J-2))*(aa-A(J-1))*(aa-A(J+1))
               A2=A2/(A2D1*A2D2*A2D3)
               A3=(aa-A(J-2))*(aa-A(J-1))*(aa-A(J))
               A3=A3/(A3D1*A3D2*A3D3)
c
               bb = A0*B(J-2) + A1*B(J-1) + A2*B(J) + A3*B(J+1)
            ENDIF 
c
            GO TO 100
         ENDIF 
      end do
  100 CONTINUE
ccc      write(2,*) 'F1, F2, F3, H1, H2, H3 =',B(J-2),B(J-1),B(J),
ccc     + A(J-2), A(J-1), A(J)
ccc      write(2,*) 'A0, A1, A2, bb =',A0,A1,A2,bb
c
      RETURN
      END 

      DOUBLE PRECISION FUNCTION 
     &         ctn2_overtone(SIGMA,PN2,PTOT,T)
C
C  Computes the absorption coefficient (in cm-1) in the 2-0 CIA band of N2
C  in air (79% N2 + 21% O2) for the wavenumber SIGMA, the temperature T, 
C  the N2 partial pressure PN2, and the total pressure PTOT
C
C   
C     SIGMA : Wavenumber in CM-1
C       PN2 : Partial pressure of N2 in ATM
C      PTOT : Total pressure of N2 in ATM
C         T : Temperature in Kelvin.
C      CTN2 : Absorption coefficient in cm-1
C
C
C NOMBRE DE VALEURS TABULEES ET PAS
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      PARAMETER (NVAL=701 , PASSIG=1.D0)
C      
      COMMON /N2DAT/SIGRF(NVAL),B0(NVAL),BETA0(NVAL)
C
      DATA T0/273.16D0/
      DATA TREF/296.D0/
C
      ctn2_overtone=0.D0
      IF ( T.GT.350.D0 ) RETURN
      IF((SIGMA.LT.SIGRF(1)).OR.(SIGMA.GT.SIGRF(NVAL)))RETURN
C 
C Interpolate within Tabulated values
      IINF=INT( (SIGMA-SIGRF(1)+0.1D-4)/PASSIG ) + 1
      IINF=MIN0(IINF,NVAL-1)
      ISUP=IINF+1
      D1ST=(1.D0/TREF)-(1.D0/T)
      BINF=B0(IINF)*DEXP(BETA0(IINF)*D1ST)
      BSUP=B0(ISUP)*DEXP(BETA0(ISUP)*D1ST)
      B=BINF+(BSUP-BINF)*(SIGMA-SIGRF(IINF))/PASSIG
C
C Multiply by N2 and Total densities in Amagat
C D(amagat)=P(atm)*(T0/T)
C
      DTOT=PTOT*(T0/T)
      DN2=PN2*(T0/T)
      ctn2_overtone=B*DN2*DTOT
      RETURN
      END

      SUBROUTINE lecn2_overtone(path_n2_cia_overtone)
C
C READs the Data for the CIA of the 2-0 band of N2 in Air
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      PARAMETER (NVAL=701)
      character path_n2_cia_overtone*128
C      
      COMMON/N2DAT/SIGRF(NVAL),B0(NVAL),BETA0(NVAL)
C      
C OUVERTURE DU FICHIER, LECTURE DE L ENTETE PUIS DES DONNEES
C
c      OPEN(UNIT=3,FILE='N2_in_air_CIA_0-2_band.dat',
      OPEN(UNIT=3,FILE=path_n2_cia_overtone,
     &          STATUS='OLD',FORM='FORMATTED')
C
      DO I=1,NVAL
         READ(3,*)SIGRF(I),B0(I),BETA0(I)
         SIGRF(I)=SIGRF(I)-5.d0
      end do
C
      CLOSE(3)
      RETURN
      END


      SUBROUTINE lecn2_fundamental(path_n2_cia_fundamental,
     &  path_h2o_n2_cia_fundamental)
C
C This routine reads the data that enable calculations
C of the N2-N2 and N2-H2O collision-induced absorptions
C in the fundamental band of N2
C
C These data have been generated as explained in the paper 
C "Indirect influence of of humidity on atmospheric emission 
C  and transmission spectra near 4 microns"
C
C Creted by J-M Hartmann, March 2018
C jmhartmann@lmd.polytechnique.fr
C
C
C Number of tabulated values
      PARAMETER (NVAL=901)
C      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      character*128 path_n2_cia_fundamental,path_h2o_n2_cia_fundamental
      COMMON/N2CIA/SIGRF(NVAL),B0air(NVAL),BETA0air(NVAL),
     &   B0H2O(NVAL),BETA0H2O(NVAL)    
C      
C Open files and read data
C
c      OPEN(UNIT=3,FILE='CT-N2.N2',
      OPEN(UNIT=3,FILE=path_n2_cia_fundamental,
     &          STATUS='OLD',FORM='FORMATTED')
C Read header then read data
      read(3,*) nlhead,ncol
      DO I=2,nlhead
         READ(3,*)
      end do
      DO I=1,NVAL
         READ(3,*)SIGRF(I),B0air(I),BETA0air(I)
      end do
      CLOSE(3)
C
c      OPEN(UNIT=3,FILE='CT-N2.H2O',
      OPEN(UNIT=3,FILE=path_h2o_n2_cia_fundamental,
     &          STATUS='OLD',FORM='FORMATTED')

C Read header then read data
      read(3,*) nlhead,ncol
      DO I=2,nlhead
         READ(3,*)
      end do
      DO I=1,NVAL
         READ(3,*)SIGRF(I),B0H2O(I),BETA0H2O(I)
      end do
      CLOSE(3)
      RETURN
      END


      DOUBLE PRECISION FUNCTION 
     &         ctn2_fundamental(SIGMA,PN2,PH2O,PTOT,T)
C
C This routine computes the absorption coefficient in the collision
C nnduced fundamental absorption band of N2 for air in the presence
C of some humidity.using the data that have been read by Subroutine "lecn2"
C
C The arguments and their units are the following
C    Sigma: wavenumber in units of "1/cm" (reciprocal centimeter)
C    PN2  : partial pressure of N2 in units of "atm" (atmosphere)
C    PH2O : partial pressure of H2O in units of "atm" (atmosphere)
C    PTOT : total pressure in units of "atm" (atmosphere)
C    T :    temperature in units of "K" (Kelvin)
C    CTN2:  absorption coefficient for the considered conditions
C           in units of "1/cm" (reciprocal centimeter). Hence, for
C           an optical path of legth L, the transmission is
C           trans = exp(-CTN2*L) where L has to be in centimeer units 
C
C Important: if the nominal N2 vmr of 0.781 is used to compute
C PN2, then PN2=0.781*(PTOT-PH2O) and NOT PN2=0.781*PTOT
C
C This routine uses a model that is described in the paper 
C "Indirect influence of of humidity on atmospheric emission 
C  and transmission spectra near 4 microns"
C
C Creted by J-M Hartmann, March 2018
C jmhartmann@lmd.polytechnique.fr
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C Number of tabulated values and wavenumber step
C      
      PARAMETER (NVAL=901 , StpSig=1.D0)
C      
C Tabulated values of the data for N2-N2 and N2-H2O
C These have been read by routine "lecn2"
C
      COMMON/N2CIA/SIGRF(NVAL),B0air(NVAL),BETA0air(NVAL),
     &   B0H2O(NVAL),BETA0H2O(NVAL)    
C
      DATA T0/273.16D0/
      DATA TREF/296.D0/
C
C
      ctn2_fundamental=0.D0
      IF ( T.GT.350.D0 ) RETURN
      IF((SIGMA.LT.SIGRF(1)).OR.(SIGMA.GT.SIGRF(NVAL)))RETURN
C 
C Compute the N2-N2 and N2-H2O CIA absorption coefficients
C (Bair and BH2o, respectively) by using the exponential Temperature
C dependence from the tabulated values and a liner interpolation versus
C wavenumber using the two sorrounding points (INF and SUP)
c  The CIA at T is computed from B0*exp[BETA0*(1/Tref-1/T)]
C
      IINF=INT( (SIGMA-SIGRF(1)+0.1D-4)/StpSig ) + 1
      IINF=MIN0(IINF,NVAL-1)
      ISUP=IINF+1
      D1ST=(1.D0/TREF)-(1.D0/T)
      BINFair=B0air(IINF)*DEXP(BETA0air(IINF)*D1ST)
      BSUPair=B0air(ISUP)*DEXP(BETA0air(ISUP)*D1ST)
      Bair=BINFair+(BSUPair-BINFair)*(SIGMA-SIGRF(IINF))/StpSig
      BINFH2O=B0H2O(IINF)*DEXP(BETA0H2O(IINF)*D1ST)
      BSUPH2O=B0H2O(ISUP)*DEXP(BETA0H2O(ISUP)*D1ST)
      BH2O=BINFH2O+(BSUPH2O-BINFH2O)*(SIGMA-SIGRF(IINF))/StpSig
C
C Then correct Bair by introducing the contribution of the N2-O2 CIA
C
      Bair=Bair*(0.79 + 0.21*(1.294D0-0.4545D0*T/TREF))
C
C Switch from pressures (in atm) to densities (in amagat)
C and compute CIA by combining dry air (N2+O2) and H2O
C contributions
C
      DTOT=PTOT*(T0/T)
      DN2=PN2*(T0/T)
      DH2O=PH2O*(T0/T)
      ctn2_fundamental=DN2*(Bair*(DTOT-DH2O)+BH2O*DH2O)
      RETURN
      END
