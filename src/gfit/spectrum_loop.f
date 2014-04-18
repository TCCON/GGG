      subroutine spectrum_loop(winfo,debug,
     & lun_col,lcl,colabel,colfile_format,lspmax,
     & runlog,akfile,rayfile,mavfile,targmol,linefiles,
     & parfile,
     & apx,apu,dplist,iptg,ipcl,ipfs,ipzo,ipcf,
     & ntg,ncbf,nfp,
     & speci,nspeci_iso,solarll,pars,sptfile)
c
c   Inputs:
c     ntg           I*4  number of target gases
c     speci(ntg)    I*4  The species # for the target gases/isotopologs
c     pars(ntg)     C**  The names of the target gases/isotopologs
c     nmp           I*4  Number of measured points (in spectrum)
c     ncbf,         I*4  number of continuum terms (basis functions)
c     nfp,          I*4  number of fited parameters = ntg+ncbf+n
c     obsrvd(nmp)   R*4  Measured spectrum (y)
c     apx(nfp)      R*4  A priori state vector (xa)
c     apu(nfp)      R*4  A priori state vector uncertainties (SQRT(diag(Sa)))
c
c   Outputs:
c     calcul(nmp)   R*4  Calculated spectrum (f(x))
c     cx(nfp)       R*4  State vector (x)
c     ex(nfp)       R*4  State vector uncertainties
c     sgshift       R*4  Solar-Gas Shift (ppm)

c  mode=0    Checks that all input files are readable before investing
c            a lot of time on doing the VAC calculation.
c  mode=1    Does full calculation

      implicit none
      include "../ggg_const_params.f"
      include "../ggg_int_params.f"
      include "const_params.f"
      include "int_params.f"

      logical  debug

      integer*4 i,lf,fbc,reclen_solarll,lr,
     & ncbf,            ! number of continuum terms (basis functions)
     & imode,           ! 0 = processing mode (do not stop)
     & lso,lsos,
     & lsn,lsnd,
     & nfov,kfov,
     & ncall,           ! counts the number of times that subroutine is called
     & oaflag,          ! Indicates whether absorbtion coefficients currently
c                         in memory are obsolete (=1) or not (=0).
     & lcolon,
     & totit,nn,ncol,nspectra,mspectra,
     & iptg,ipcl,ipfs,ipzo,ipcf,
     & mcp,ncp,
     & mva,nva,         ! The number of precomputed anscoption coefficients
c     & mvmr,
     & istat,ifm,
     & nscycle,
     & mspt,ispec,
     & freq_flag,       ! =1  presents spectral fits in solar rest frame. 
                        ! =0  presents spectral fits in atmosphere rest frame.
     & lun_ak,
     & lun_sts,         ! Solar Transmittance Spectrum
     & lnbc,
     & nspeci_iso,jspeci,
     & nhw,ldec,nmpfp,
     & nmp,imp,
     & mii, nii,            ! Dimension of interpolated (by factor LDEC) ILS
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
     & mit,nit,
     & iyr,iset,
     & kspflag,ifirst,ilast,bytepw,possp,
     & kcp1,kcp2,
     & nsh,nhwmax,
     & lun_col,lun_spt,lunr_mav,lunr_ray,lun_rlg,mav_count

      parameter (mva=166000000,mcp=1650000,
     & nscycle=25,
c     mvmr=28000,
     & mii=103847,
     & mslpd=10*mmp,
     & mspxv=12*mcp)

      parameter (lun_sts=24,lun_rlg=25,lun_ak=26,
     & lunr_ray=27,lunr_mav=28,lun_spt=29)

      integer*4
     & speci(ntg),
     & targmol(nspeci_iso)

      real*4
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
     & cont_level,cont_tilt,cont_curv,xzo,
     & dd,
     & sza_ray,sza_raywas,bend,
     & tot,
     & vmr(mspeci*mlev),
     & vpf(mspeci*mlev),
     & spver(mlev),splos(mlev),
     & cp(mlev),
     & z(mlev),t(mlev),p(mlev),d(mlev),
     & corrld,           ! Factor=[apodixed resolution]/[spectral point spacing
     & solar_gas_shift,sgshift,
c     & fovcf,
     & tbar, rnoise,detnoise,sphnoise,
     & slpd(mslpd),
     & pd((mmp+mfp)*mfp), spts(mcp),
     & spxv(mspxv),
     & solzen,roc,fbar,
     & xfs,rdum

      real*8
     & dopp,
     & riair,
     & tottc,tottc2,toterr,avgtc,rmstc,avgcl,avgrms,
     & sssss,            ! Spectral offset expressed as fraction of GINT.
     & frac,             ! fractional size of FOVO compared with solar diameter
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
     & fzero,            ! frequency of the zero'th computed VAC
     & startm,           ! frequency of first returned measured point
     & nus,              ! microwindow starting frequency (cm-1)
     & nue,              ! microwindow ending frequency (cm-1)
     & graw,             ! spacing of raw spectrum (cm-1) from GETINFO
     & gint,             ! spacing of OBSRVD (cm-1) after interpolation
     & grid,             ! spacing of primitive spectrum 
     & vbar,             ! mean frequency of measured spectrum
     & hwid              ! half-width of measured spectrum

      real*8 
     & oblat,            ! observation latitude (deg).
     & oblon,            ! observation longitude (deg).
     & obalt,            ! observation altitude (km)
     & zpdtim,           ! Time of ZPD (UT hours)
     & asza,             ! astronomical solar zenith angle (unrefracted)
     & zenoff,           ! zenith pointing offset
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

      character winfo*(*),specname*(nchar),pars(ntg)*(*),
     & data_fmt_read_rl*256,
     & col_labels_rl*320,
     & sptfile*(*),akpath*(mfilepath),akfile*(*),specpath*(mfilepath),
     & sptpath*(mfilepath),header*256,rayfile_format*256,
     & solarll*(*),runlabmav*(nchar),
     & oformat*12,colabel*(*),
     & mavstring*64,linefiles*(*),parfile*(*),dplist*(*),
     & col1*1,apf*2,rayfile*(*),specray*(nchar),runlog*(*),mavfile*(*),
c     & ss(nfp)*4,
     & string*80,colfile_format*(*)

      parameter (mspectra=999999,resmax=0.375d0,opdmin=0.5d0/resmax)
c

      save ispec,ncall
      data ispec/0/
      data ncall/0/

      rayfile_format='(a,200f12.0)'

      oaflag=1
      mavflag=0
      imode=0
c  Read max # of SPT files (if a value is provided on the SPT line of the .ggg file)
      mspt=3500  ! default value
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
         write(lun_col,oformat)' Spectrum            ',
     &   colabel(:lnbc(colabel))
      endif 

      open(lunr_mav,file=mavfile,status='old')
      read(lunr_mav,*)
      read(lunr_mav,'(14x,a)')string
      lcolon=index(string,':')   ! Next Spectrum:
      read(string(lcolon+1:),'(a)') runlabmav   ! First Next Spectrum

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
 
c      if(nlev_ray*nspeci_iso.gt.mvmr) then
c         write(*,*)'nlev_ray*nspeci_iso,mvmr=',nlev_ray*nspeci_iso,mvmr
c         stop ' spectrum_loop: Increase parameter mvmr'
c      endif
c
      read(winfo,*) frqcen,width,mit,defapo,interp,freq_flag
      grid=0.666666d-06*frqcen
      nus=frqcen-width/2
      nue=frqcen+width/2
12    kcp1=int(nus/grid)
      kcp2=int(nue/grid)
      nsh=int(2+2*resmax/grid)
      nhwmax=nint(nscycle*resmax/grid)
      ncp=kcp2-kcp1+2*nhwmax+2*nsh
      fzero=grid*(kcp1-nsh-nhwmax)
      nva=ncp*nlev_ray*(ntg+1)
      if(nva+ncp.gt.mva) then
        write(6,*)'Increase MVA from',mva,' to',nva+ncp
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
         stop 'next_spectrum: Increase parameter MSPXV'
      endif
C
      avgrms=0.0d0
      tottc=0.0d0
      tottc2=0.0d0
      toterr=0.0d0
      avgcl=0.0d0
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
      if(debug) write(*,*)' Main loop...',nspectra,mspectra
      do ispec=1,mspectra         !  Main fitting loop over spectra
141     call read_runlog_data_record(lun_rlg,data_fmt_read_rl,
     &  col1,specname,iyr,iset,zpdtim,
     &  oblat,oblon,obalt,asza,zenoff,azim,osds,
     &  opd,fovi,fovo,amal,ifirst,ilast,graw,possp,bytepw,zoff,snr,apf,
     &  tins,pins,hins,tout,pout,hout,
     &  sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)
        if(debug) write(*,*) 'runlab, istat=',specname, istat
        if(istat.ne.0) exit
        if(col1.eq.':') cycle

        opdmax=dabs(0.5/graw)
        lsp=lnbc(specname)
        specname=specname(:lsp)
c
        read(lunr_ray,rayfile_format) specray,rdum,rdum,sza_ray,bend,
     &  rdum,zmin,(splos(j),j=1,nlev_ray)
c        write(*,*) 'zmin=',zmin
c        write(*,*)'nlev_ray,splos=',nlev_ray,specray,(splos(j),j=1,nlev_ray)
c        write(37,*)zmin,zminwas,sza_ray,sza_raywas,
c     &  (zmin-zminwas)/(sza_ray-sza_raywas)
        zminwas=zmin
        sza_raywas=sza_ray
        specray=specray(:lsp)
        if(specray.ne.specname) then
           write(6,'(a,2x,a)') specray,specname
           stop 'spectrum mismatch 2'
        endif

c  If MIT > 0, check that requested spectrum is on disk,
c  and that it covers the specified spectral interval.
        kspflag=0
        call gindfile(dplist,specname,specpath)
        if(lnbc(specpath).eq.0) kspflag=2
        if(debug) write(*,*)' kspflag = ',
     &     kspflag,'  '//specname,specpath
c
c  Apply air-to-vacuum correction
c        write(*,*)'tins,pins,hins',tins,pins,hins
        if(kspflag.lt.2) graw=graw*riair(lasf,tins,pins,hins)/
     &  riair(frqcen,tins,pins,hins)

c  Apply FOV corrections
c        graw=graw*(1.D0+(amal**2+fovi**2)/16)  ! FOV correction

        resn=0.5d0/opd
        if(opd.gt.opdmax) then
           resn=0.5d0/opdmax
           write(*,*)'Warning: OPD out of range: reset to OPDmax: ',
     &     specname(:lsp),opd,opdmax
           if (imode.gt.0) stop 'OPD > OPDmax'
        endif
        if(opd.lt.opdmin) then
           resn=0.5d0/opdmin
           write(*,*)'Warning: OPD out of range: reset to OPDmin: ',
     &     specname(:lsp),opd,opdmin
           if (imode.gt.0) stop 'OPD < OPDmin'
        endif

c        if(1.0001*resn.lt.graw) then
c           write(*,*) 'resn,graw=',resn,graw,specname
c           stop ' resn < graw'
c        endif

c  Measured spectrum must be wider than fitting interval to
c  allow convolution with apodizing/interpolating ILS
        dd=nscycle*resn ! half-width of the slit function in cm-1 (always +ve)
        vbar=0.5d0*graw*(ilast+ifirst)
        hwid=0.5d0*dabs(graw*(ilast-ifirst))
c        write(*,*)ilast,ifirst,vbar,hwid,dd
        if(kspflag.eq.2) dd=0.0                 ! Making synthetic spectrum
        if(debug)write(*,*)kspflag,nus,nue,vbar,hwid,vbar-hwid,vbar+hwid
        if( nus-dd .lt. vbar-hwid ) kspflag=1   ! Lower window limit < disk file
        if( nue+dd .gt. vbar+hwid ) kspflag=1   ! Upper window limit > disk file
c========================================================================
c  Read model & vmr information (SUNRUN.MAV)
c      lsp=lnbc(specname)
        if(kspflag.eq.0 .and. lsp.gt.lspmax) lspmax=lsp  ! Longest fitted spectrum name
c      write(*,*)'specname, runlabmav = ',
c     & specname(:24),runlabmav(:24),kspflag,oaflag
      if(specname.eq.runlabmav) then
        call read_mavfile_body(lunr_mav,nlev_ray,nspeci_iso,
     &  z,t,p,d,vmr,parfile)
         mavflag=1
         read(lunr_mav,'(a)',end=66) mavstring
         if(mavstring(1:14).eq.'Next Spectrum:') then
            read(mavstring(15:),'(a)') runlabmav
         else
            write(*,*) mavstring
            write(*,*)'Failed to find Next Spectrum: key-word'
            stop
         endif
66       continue
         if(index(winfo,' sa_temp ').gt.0)
     &   call vadd(t,1,5.,0,t,1,nlev_ray) 
         if(index(winfo,' sa_pres ').gt.0)
     &   call vmul(p,1,.95,0,p,1,nlev_ray)
         oaflag=1  ! Absorption coefficients currently in memory are obsolete
      endif  ! specname.eq.runlabmav

      if(mavflag.le.0) then
      write(*,*)'Mismatch between spectrum names in mavfile and runlog'
      write(*,*) runlabmav,specname
      stop 'mavfile data have not been read in -- cannot proceed.'
      endif

      if(kspflag.eq.0 .and. oaflag.ge.1) then  !  current spectrum encompasses window
c  Pre-compute absorption coefficient
c         nva=ncp*nlev_ray*(ntg+1)
c         write(*,*)' Calling abscoj...', kspflag
         if(ncall.ge.1) then  ! First call, skip time-consuming stuff
         call vmov(zero,0,vac,1,nva)
         call abscoj(specname,nlev_ray,t,p,d,nspeci_iso,targmol,vmr,vpf,
     &   linefiles,parfile,fzero,grid,ncp,vac,vac(nva+1))
         write(6,oformat)'  Spectrum      ',colabel(:lcl+48)
         endif          ! (ncall.ge.1) then
c         write(*,*)' Called abscoj...'
         mav_count=mav_count+1
         oaflag=0   ! Absorption coefficient in memory are up to date
      endif          ! (kspecflag.eq.0 .and. oaflag.eq.1) then 


      if((mit.gt.0.and.kspflag.gt.0).or.(mit.eq.0.and.kspflag.eq.1))then
         if(debug) then
         write(*,*)' Requested interval: ',nus,' to ',nue
         write(*,*)' Required interval: ',nus-dd,' to ',nue+dd
         write(*,'(a,a,f9.3,a4,f9.3,a5)') specname,
     & ' only encompasses',ifirst*graw,' to ',ilast*graw,' cm-1'
         endif
         goto 141 ! skip missing/partial spectrum
      endif
c=========================================================
c      resn=0.5d0/opd
c      if(resn.lt.graw) resn=graw
      if(index(winfo,' sa_fovi ').gt.0) fovi=fovi*1.07
      rect=frqcen*(fovi**2+amal**2)/8  ! old code
      if(ncall.ge.1) then
c
c  This code from DG is a crude attempt to simulate
c  a frequency-independent misalignment effect.
c      if(runlog(ir+1:ir+2).eq.'in') then
c         rect=frqcen*(fovi**2)/8 + 2100.D0*(amal**2)/8
c      elseif(runlog(ir+1:ir+2).eq.'hg') then
c         rect=frqcen*(fovi**2)/8 + 800.D0*(amal**2)/8
c      else
c         rect=frqcen*(fovi**2+amal**2)/8
c      endif
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
     &       specname,apf,' ???  Unknown apodization function'
             stop
           endif
        endif

c---------------------------------------------------------
c  FIND the spectrum file, return the PATH to the spectrum
      call vmov(zero,0,obsrvd,1,mmp)
c      write(*,*)'kspflag=',kspflag
      if(kspflag.eq.0) then
        call jetspe(specpath,resn,graw,ifirst,ilast,possp,bytepw,nus,
     &  nue,apo_m,interp,zero,zero,
     &  obsrvd,mmp,nmp,startm,gint,rc)
        if(debug) write(*,*)'obsrvd=',obsrvd(1),obsrvd(2),obsrvd(3)
        if(rc.ne.0) then
           write(6,*)' Error in JETSPE. Spectrum ',specname,rc,nmp
           write(6,*)' This error should never happen'
           kspflag=1
        endif
      endif  ! kspflag.eq.0
      if(kspflag.eq.2) then
         startm=nus
         gint=graw/interp
         nmp=1+(nue-nus)/gint
      endif
      if(nmp.gt.mmp) then
        write(*,*)'nmp,mmp=',nmp,mmp
        stop 'Increase parameter MMP'
      endif

      lsn=index(winfo,' sa_snr ')
      lsnd=index(winfo,' sa_snr/')
      if(lsn.gt.0 .or. lsnd.gt.0) then  ! Add noise
c  Add noise and systematic error to OBSRVD (assumes a continuum level of 1.0)
      call vdot(obsrvd,1,unity,0,tot,nmp)
      tbar=tot/nmp
      
      detnoise=0.0010 ! MATMOS detector noise
      sphnoise=0.0025 ! MATMOS source photon noise
c      detnoise=0.0006 ! MATMOS detector noise
c      sphnoise=0.0015 ! MATMOS source photon noise

      if(lsnd.gt.0) then
        read(winfo(lsnd+8:),*) fsn ! =1.5 for NOMAD and ACS (QQQQQ)
        detnoise=detnoise/fsn  !  detnoise=0.0004 ! NOMAD detector noise
        sphnoise=sphnoise/fsn  !  sphnoise=0.0010 ! NOMAD source photon noise
      endif

c  Add detector noise and source photon noise in quadrature
      rnoise=SQRT(tbar*sphnoise**2+detnoise**2)
c      write(*,*)'tbar,rnoise,SNR =',tbar,rnoise,tbar/rnoise
      do imp=1,nmp
        obsrvd(imp)=obsrvd(imp)+rnoise*gasdev(iseed) ! Add random noise
c      obsrvd(imp)=obsrvd(imp)+0.00012*gasdev(iseed)*opd ! add 0.4% rms random noise (iATMOS SNR=250)
c      obsrvd(imp)=obsrvd(imp)+0.0002*gasdev(iseed)*opd ! add 0.5% rms random noise
c      obsrvd(imp)=obsrvd(imp)+0.0001*gasdev(iseed)*opd ! add 0.25% rms random noise (MATMOS SNR=400)
c      obsrvd(imp)=obsrvd(imp)*(1+0.01*4*imp*(nmp-imp)/nmp/nmp)  ! add +- 1% Cont Curv
c      obsrvd(imp)=obsrvd(imp)*(1+0.005*obsrvd(imp))  !  0.5%  detector non-linearity
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
       obsrvd(imp)=obsrvd(imp)+0.001*tbar*cos(2*dpi*imp*gint/0.41)
       obsrvd(imp)=obsrvd(imp)+0.001*tbar*cos(2*dpi*imp*gint/1.01)
      end do
      endif  !  Add channelling

      nmpfp=nmp+nfp
c  Pre-compute ILS (oversampled by a factor LDEC) and normalize each
c  interleaved component to unity.
      rdec=gint/grid
      nhw=nint(nscycle*resn/grid)
      ldec=1+int(16*grid/resn)
      nii=1+2*nhw*ldec
      resnog=ldec*resn/grid
      rectog=ldec*rect/grid
      call profzl(apo_c,nii,resnog,rectog,0.0d0,slit)
c      write(*,*)'apo_c=',apo_c,apo_m,resn, rect, ldec, grid
c      write(*,*)resnog, rectog, nhw, nii
      do j=ldec,1,-1
        call vdot(slit(j),ldec,unity,0,tot,2*nhw)
        call vmul(slit(j),ldec,1.0/tot,0,slit(j),ldec,2*nhw)
      end do
      slit(nii)=slit(nii)/tot

c Write out ILS
c      open(56, file='fort.56')
c      write(56,*) 2,2
c      write(56,*)' v s'
c      do j=1,nii,ldec
c         write(56,*)grid*(j-float(nii+1)/2)/ldec,slit(j)
c      end do
c      close(56)
c===================================================================
c  Error estimation naively assumes spectral points are linearly independent.
c  CORRLD is an estimate of how dependent neighbouring points really are.
c  CORRLD = 1 means that points are truly linearly independent.
c  Variances must be scaled by CORRLD to ensure that when the same spectrum is
c  analyzed under different values of APO or INTERP, the same size errors result
c      corrld=sqrt(rect**2+((1.+0.4*apo_c)*resn)**2+(resn-0.5d0/opd)**2)
c     & /gint
      corrld=sqrt(rect**2+((1.+0.4*apo_c)*resn)**2)/gint
c
c  SSSSS is the frequency offset expressed as a fraction
c  of the measured spectral spacing.
      sssss=sngl((startm/grid-kcp1+nsh-1+nhwmax-nhw)/rdec)

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
         do ilev=1,nlev_ray
            call vsma(vac(jva),1,ckm2cm*splos(ilev),spxv(jsp),1,
     &      spxv(jsp),1,ncp)
c            write(*,*)'spxv:',ilev,z(ilev),jva,splos(ilev),vac(jva)
            jva=jva+ncp
         end do
c         if(jtg.eq.1) write(44,'(3f8.3,1p7e10.2)')zmin,sza_ray,bend,
c     &   (-spxv(jsp+jj),jj=0,ncp-1,ncp/6-1)
         jsp=jsp+ncp
      end do

c   Add solar optical thickness spectrum to non-target ones in SPXV(JSP)
c   Use SPXV(JSP+NCP) as work space for Voigt functions.
      lso=index(winfo,' so ')
      lsos=index(winfo,' so/')
      if( lso .gt. 0 .or. lsos .gt. 0) then
         frac=fovo/9.2e-3  ! the sun is 9.2 mrad in diameter on average
         if(index(winfo,' sa_osds ').gt.0) osds=osds-1.0  ! Add 1ppm
         lr=lnbc(solarll)
         read(solarll(lr-2:lr),*) reclen_solarll
         if(reclen_solarll.eq.108) then
           call solar_pts(lun_sts,solarll,
     &     fzero*(1.d0+osds*1.E-6),grid*(1.d0+osds*1.E-6),frac,spts,ncp)
         else
           stop 'Upgrade solar linelist to .108 format'
         endif
         if(lsos.gt.0) then
           read(winfo(lsos+4:),*) ffr  ! absorption reduction factor
           write(*,*)'Reducing solar absorption depths by factor',ffr
           do i=1,ncp
             spts(i)=1.0-(1.0-spts(i))/ffr
           end do
         endif
      else
         osds=0.0d0
         call vmov(unity,0,spts,1,ncp)
      endif  ! lso .gt. 0
      if(index(winfo,' sa_zoff ').gt.0) zoff=zoff+0.005

c===================================================================
c  Compute the vertical slant path distances
c      write(*,*)'nlev_ray,ncell=',nlev_ray,ncell
      if (nlev_ray.le.ncell+1) then
         spver(1)=splos(1)
         spver(nlev_ray)=splos(nlev_ray)
      else
c         write(*,*)'Valling CVP'
         call compute_vertical_paths(ncell,zmin,z,d,spver,nlev_ray)
      endif  ! nlev_ray.le.2
c      do j=1,nlev_ray
c      write(*,*) j,z(j),d(j),spver(j),splos(j)
c      end do
c      write(*,*)'splos=',(splos(j),j=1,nlev_ray)
c      write(*,*)'spves=',(spver(j),j=1,nlev_ray)
c===================================================================
c Compute vertical and LOS slant columns for the target gases
c Do in reverse order so that the concentration profile of the
c first target gas is retained in CP at the end.
      do jtg=ntg,1,-1
         jspeci=speci(jtg)
         call vmul(vmr(jspeci),nspeci_iso,d,1,cp,1,nlev_ray)
         call vdot(cp,1,spver,1,overcol(jtg),nlev_ray)
         call vdot(cp,1,splos,1,oloscol(jtg),nlev_ray)
         overcol(jtg)=overcol(jtg)*ckm2cm
         oloscol(jtg)=oloscol(jtg)*ckm2cm
c      write(*,*)'jtg,jspeci,overcol,oloscol',jtg,jspeci,
c     & overcol(jtg),oloscol(jtg)
      end do
c      do i=1,nlev_ray
c        write(*,*)'vmr=',i,z(i),(vmr(speci(jtg)+nspeci_iso*(i-1)),jtg=1,ntg)
c      end do
c
      if(ntg.eq.0) then
c   compute total column amounts for air (vmr=1.000)
         call vdot(d,1,spver,1,overcol(1),nlev_ray)
         call vdot(d,1,splos,1,oloscol(1),nlev_ray)
c      write(*,*)'ntg,overcol,oloscol',ntg,overcol(1),oloscol(1)
      endif
c=====================================================================
c  The structure of the state vector is as follows. Elements:
c      1 to NTG            are the Target Gases
c  NTG+1 to NTG+NCBF       are the Continuum Basis Funcbfions Coefficients
c  NTG+NCBF+NFS            is the FS (if NFS>0)
c  NTG+NCBF+NFS+NZO        is the ZO (if NZO>0)
c  
c  Read the initial values of the parameters to be fitted and their estimated
c  a priori variances.
c      open(lun_apx,file=ap_file,status='old')
c      read(lun_apx,*)
c      kcbf=ntg+1
c      do j=1,5       !  Read a priori values & uncertainties
c         read(lun_apx,'(f5.0,f9.0,1x,a4)') apx(kcbf),apu(kcbf),ss(kcbf)
c         if(index(winfo,ss(kcbf)).gt.0) kcbf=kcbf+1
c      end do
c      do jj=ntg+1,ntg+5
c         read(lun_apx,'(f5.0,f9.0,1x,a4)') apx(jj),apu(jj),ss(jj) 
c             if(ss(jj).eq.' cl ') then
c             ipcl=jj
c         elseif(ss(jj).eq.' ct ') then
c             ipct=jj
c         elseif(ss(jj).eq.' cc ') then
c             ipcc=jj
c         elseif(ss(jj).eq.' fs ') then
c             ipfs=jj
c         elseif(ss(jj).eq.' zo ') then
c             ipzo=jj
c         else
c             write(*,*) ntg,jj,jj-ntg,ss(jj)
c             stop 'Unknown/Unreadable mnemonic found in a priori file'
c         endif
c      end do

c      read(lun_apx,*) apx(ipcl),apu(ipcl),ss(ipcl) ! Continuum Level
c      read(lun_apx,*) apx(ipct),apu(ipct),ss(ipct) ! Continuum Tilt
c      read(lun_apx,*) apx(ipcc),apu(ipcc),ss(ipcc) ! Continuum Curvature
c      read(lun_apx,*) apx(ipfs),apu(ipfs),ss(ipfs) ! Frequency Shift
c      read(lun_apx,*) apx(ipzo),apu(ipzo),ss(ipzo) ! Zero Offset
c      read(lun_apx,*)                              ! Solar Scaling
c      if(ntg.ge.1) read(lun_apx,*) apx(1),apu(1)   ! First Target Gas
c      if(ntg.ge.2) read(lun_apx,*) apx(2),apu(2)   ! Other Target Gases
c      close(lun_apx)
c      call vmov(apx(2),0,apx(3),1,ntg-2) ! absorber amount
c      apx(ipzo)=apx(ipzo)+sngl(zoff)
c      call vmov(apu(2),0,apu(3),1,ntg-2)  ! vmr of non-target gases
c      apu(ipfs)=apu(ipfs)*sngl(1.d0+rdec) ! shift   up from 0.5 4-DEC-98 gct
c==========================================================================
c      write(*,*)' ntg=',ntg
c      write(*,*)'sl: apx(ipct)=',ipct,apx(ipct)
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
      call vmov(aprx,1,cx,1,nfp)  ! Initialize to A PRIORI uncertainties each spectrum
      if(debug)
     &  write(6,*)  'It    CL       CT    CC    FS    ZOFF  RMS/CL'//
     &'  Vfact1  Vfact2  Vfact3  Vfact4  Vfact5  Vfact6'

      if(index(winfo,' prof_ret_only ').gt.0) ifm=4  ! Profile Retrieval only
      if(index(winfo,' prof_ret+scale ').gt.0) ifm=5  ! Profile Retrieval
      kfov=index(winfo,'nfov=')
      if(kfov.gt.0) then
         read(winfo(kfov+5:),'(i1)')nfov
         ifm=3
         solzen=sngl(asza+zenoff)
         roc=6378.
         fbar=(ifirst+ilast)*graw/2
c         write(*,*) 'Calling do_retrieval3...'

         call do_retrieval3(obsrvd,nmp,aprx,apru,slit,nii,
     &   iptg,ipcl,ipfs,ipzo,ipcf,xzo,
     &   cont_level,cont_tilt,cont_curv,
     &   nfov,z,t,p,solzen,fovo,roc,obalt,wavtkr,fbar,
     &   ldec,rdec,spts,spxv,vac,splos,nlev_ray,ncp,ntg,ncbf,nfp,snr,
     &   corrld,sssss,debug,mit,nit,cont,calcul,rmsocl,
     &   cx,ex,slpd,pd,ssnmp)
      else  !  FOV center only.
         ifm=1
         nfov=1
         if(debug)write(*,*)'Calling do_retrieval...'
         call do_retrieval(obsrvd,nmp,aprx,apru,slit,nii,
     &   iptg,ipcl,ipfs,ipzo,ipcf,xzo,
     & cont_level,cont_tilt,cont_curv,
     &   ldec,rdec,spts,spxv,vac,splos,nlev_ray,ncp,ntg,ncbf,nfp,snr,
     &   corrld,sssss,debug,mit,nit,cont,calcul,rmsocl,cx,ex,
     &   slpd,pd,ssnmp)
         if(debug)write(*,*)'Called do_retrieval.'
      endif

c      if(ifm.eq.1) then   ! FOV center ray only
c      elseif(ifm.eq.3) then  ! Full numerical FOV integration (slow)
c      endif

c  Set error bars very small for air, otherwise it will
c  dominate the CO2/Air ratio uncertainties.
      if(index(winfo,' air ').gt.0) ex(1)=1.0e-08

c  Write .spt file for the first MSPT spectral fits
c       write(*,*)'sptpath=',sptpath, ispec,mspt,100*rmsocl
      if(ipfs.gt.0) then
         xfs=cx(ipfs)
      else
         xfs=0.0
      endif
      if(freq_flag.eq.0) then
         dopp=xfs*gint/((nus+nue)/2)
      else
         dopp=osds*1.E-06
      endif
      if(ispec .le. mspt) then
         sptpath=sptfile(:lf-1)//specname
         call write_spt(lun_spt,winfo,sptpath,
     &   obsrvd,calcul,cont,cx,ex,startm+gint*(xfs),
     &   dopp,gint,overcol,pars,asza+zenoff,obalt,zmin,
     &   xzo,abs(100*rmsocl),frac,pd,ssnmp,nmp,nmpfp,ntg,nfp)
      endif

c  Solar_Gas_shift is the difference of the solar shift and the gas shift
c  expressed in terms of the observed spectral point spacing
c  Multiply by GINT to convert to cm-1.
c   Divide by frequency to normalize to a dimensionless stretch which
c  should then  be the same for all windows. Multiply by 10^6 for ppm.
      if(index(winfo,' so ').gt.0) then
         sgshift=1E6*2*gint/frqcen*
     &   solar_gas_shift(cont_level,cont_tilt,obsrvd,calcul,ssnmp,nmp)
      else
         sgshift=0.0
      endif

c  Write the state vector and uncertainties to the .col file
      call write_col(lun_col,colfile_format,specname,lspmax,nit,
     & ncbf,xzo,xfs,cont_level,cont_tilt,cont_curv,
     & rmsocl,sgshift,gint,zmin,oloscol,overcol,cx,ex,ntg,nfp)
c  And to the screen
      call write_col(6,colfile_format,specname,lspmax,nit,
     & ncbf,xzo,xfs,cont_level,cont_tilt,cont_curv,
     & rmsocl,sgshift,gint,zmin,oloscol,overcol,cx,ex,ntg,nfp)

c      if(debug .and. nit.eq.mit+1) stop 'nit=mit+1' 

c  Output PD's (weighting functions), for subsequent use
c  in deriving averaging kernels.
         if( index(winfo,' ak ') .gt. 0) then
            akpath=akfile(:lnbc(akfile))//specname
            open(lun_ak,file=akpath,status='unknown')
            call fm(lun_ak,slit,nii,
     &      iptg,ipcl,ipfs,ipzo,ipcf,
     &      ldec,spts,spxv,
     &      vac,splos,nlev_ray,ncp,rdec,sssss,
     &      cont_level,cont_tilt,cont_curv,xzo,
     &      cx,ntg,ncbf,nfp,cont,
     &      calcul,slpd,pd,nmp)
            ynoise=2.5*cont_level*corrld/sngl(0.1d0+snr)
c  Skip levels representing the cells, so start ilev at ncell+1.
            write(lun_ak,*)(splos(ilev)*cp(ilev),ilev=ncell+1,nlev_ray)
            write(lun_ak,*) pout/1013.25
            write(lun_ak,*)(p(ilev),ilev=ncell+1,nlev_ray)
            write(lun_ak,*)(ynoise/apru(ifp),ifp=1,nfp)
            write(lun_ak,*)((aprx(ifp)-cx(ifp))*ynoise/apru(ifp),
     &      ifp=1,nfp)
            close(lun_ak)
         endif

      totit=totit+nit
      avgrms=avgrms+1
      avgcl=avgcl+abs(1/rmsocl)
      toterr=toterr+1.0/ex(1)**2
      tottc=tottc+cx(1)/ex(1)**2
      tottc2=tottc2+((cx(1)-1)/ex(1))**2

      endif        !  (ncall.ge.1) then
      end do   !  ispec=1,mspectra     Main fitting loop over spectra
      nspectra=ispec-1
      if(ncall.gt.0) then
         write(*,*)' Grid=',grid,'cm-1' 
         write(*,*)' Used ',float(nva+ncp)/mva,' of allocated memory'
         write(6,*)' Total number of iterations =',totit
         write(6,*)' Total number of spectra fitted  =',nspectra
         write(6,*)' Total number of mavblocks found =',mav_count
         write(6,*)' Average % RMS fit =',100*avgrms/avgcl
         avgtc=tottc/toterr
         rmstc=sqrt(abs(tottc2/toterr-(avgtc-1.)**2))
         write(6,*)' Average VSF =',avgtc,' +- ',rmstc
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
