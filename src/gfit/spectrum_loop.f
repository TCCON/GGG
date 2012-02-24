      subroutine spectrum_loop(ap_file,winfo,
     & lun_col,lcl,colabel,colfile_format,
     & runlog,akfile,rayfile,mavfile,targmol,linefiles,
     & parfile,dplist,ntg,speci,nspexi,solarll,pars,sptfile)
c
c   Inputs:
c     nmp           I*4  Number of measured points (in spectrum)
c     ntg           I*4  number of target gases
c     nfp           I*4  number of fitted parameters (in state vector)
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

      logical
     & debug

      integer*4 i,nlhead,lf,fbc,reclen_solarll,
     & lso,lsos,
     & lsn,lsnd,
     & ncall,           ! counts the number of times that subroutine is called
     & lcolon,
     & totit,nn,ncol,nspectra,mspectra,
     & n1,n2,n3,n4,n5,
     & mcp,ncp,
     & mva,nva,         ! The number of precomputed anscoption coefficients
     & mvmr,
     & istat,ifm,
     & nscycle,
     & mspt,ispec,
     & freq_flag,       ! =1  presents spectral fits in solar rest frame. 
                        ! =0  presents spectral fits in atmosphere rest frame.
     & lun_ak,
     & lun_sts,         ! Solar Transmittance Spectrum
     & lun_apx,         ! A priori state vector file
     & lnbc,
     & nspexi,jspeci,
     & nmp,imp,
     & nhw,ldec,nmpfp,
     & mii, nii,            ! Dimension of interpolated (by factor LDEC) ILS
     & ntg,jtg,
     & nfp,ifp,
     & jsp,jva,
     & lcl,lc,lr,lrmax,
     & interp,
     & defapo,apo_c,apo_m,rc,j,
     & mspxv,
c     & mslpd,
     & iseed,
     & nlev,ilev,
     & mit,nit,
     & iyr,iset,
     & kspflag,ifirst,ilast,bytepw,possp,
     & kcp1,kcp2,
     & nsh,nhwmax,
     & lun_col,lun_spt,lun_mav,lun_ray,lun_rlg,mavfound

      parameter (mva=160000000,mcp=1440000,
     & nscycle=25,mvmr=28000,mii=103847,
c     & mslpd=10*mmp*mtg,
     & mspxv=5*mcp)

      parameter (lun_apx=23,lun_sts=24,lun_rlg=25,lun_ak=26,
     & lun_ray=27,lun_mav=28,lun_spt=29)

      integer*4
     & speci(mtg),
     & targmol(nspexi)

      real*4
     & vac(mva),
     & ssnmp(mmp),
     & slit(mii),cx(ntg+5),ex(ntg+5),
     & obsrvd(mmp),calcul(mmp),
     & overcol(mtg),oloscol(mtg),
     & zmin,zminwas,
     & fsn,ffr,
     & gasdev,  ! random number generator
     & apx(ntg+5),apu(ntg+5),
     & ynoise,
     & rms, 
     & dd,
     & sza_ray,sza_raywas,bend,
     & tot,
     & vmr(mvmr),
     & vpf(mvmr),
     & spver(mlev),splos(mlev),
     & cp(mlev),
     & z(mlev),t(mlev),p(mlev),d(mlev),
     & corrld,           ! Factor=[apodixed resolution]/[spectral point spacing
     & solar_gas_shift,sgshift,
c     & fovcf,
     & tbar, rnoise,detnoise,sphnoise,
c     & slpd(mslpd),
     & pd((mmp+mfp)*mfp), spts(mcp),
     & spxv(mspxv),
     & solzen,roc,fbar,
     & rdum

      real*8
     & riair,
     & tottc,tottc2,toterr,avgtc,rmstc,avgcl,avgrms,
     & sssss,
     & frac,             ! fractional size of FOVO compared with solar diameter
     & resn,             ! 0.5d0/opd = half width of SINC function in cm-1
     & resmax,           ! maximum allowed value of resn
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

      character winfo*(*),ap_file*(*),specname*(nchar),pars(ntg)*(*),
     & sptfile*(*),akpath*(mfilepath),akfile*(*),specpath*(mfilepath),
     & sptpath*(mfilepath),
     & solarll*(*),runlabmav*(nchar),oformat*10,colabel*(*),
     & mavstring*64,linefiles*(*),parfile*(*),dplist*(*),
     & col1*1,apf*2,rayfile*(*),specray*(nchar),runlog*(*),mavfile*(*),
     & string*80,colfile_format*(*)

      parameter (mspectra=999999,resmax=0.375d0)
c
      save lrmax,ispec,ncall
      data lrmax/0/
      data ispec/0/
      data ncall/0/

      n1=ntg+1
      n2=ntg+2
      n3=ntg+3
      n4=ntg+4
      n5=ntg+5
      nfp=ntg+5

c  Read max # of SPT files (if a value is provided on the SPT line of the .ggg file)
      mspt=1200  ! default value
      lf=fbc(sptfile)
      if(lnbc(sptfile).gt.lf) read(sptfile(lf:),*) mspt

      lc=index(winfo,':')
      call substr(winfo(lc+1:),pars,mtg,ntg)
      if(ntg.gt.mtg) then
          write(*,*)' spectrum_loop: Error: NTG > MTG ',ntg,mtg
          stop 'Increase parameter MTG inspectrum_loop.f '
      endif

      if( index(winfo,'debug') .gt. 0 ) then
         debug=.true.
      else
         debug=.false.
      endif

      if (debug) write(*,*) 'Entered spectrum_loop: ncall=',ncall

      open(lun_mav,file=mavfile,status='old')
      read(lun_mav,*)
      read(lun_mav,'(14x,a)')string
      lcolon=index(string,':')
      read(string(lcolon+1:),'(a)') runlabmav

      open(lun_ray,file=rayfile,status='old')
      read(lun_ray,*)nn,ncol
      call skiprec(lun_ray,nn-1)
      nlev=ncol-7  ! First 7 columns are (spec,Zobs,Pobs,SZA,Bend,FOV,Zmin)
c
      if(nlev.gt.mlev) then
         write(*,*) 'nlev,mlev=',nlev,mlev
         stop 'Increase parameter MLEV'
      endif
 
      if(nlev*nspexi.gt.mvmr) then
         write(*,*)'nlev*nspexi,mvmr=',nlev*nspexi,mvmr
         stop ' spectrum_loop: Increase parameter mvmr'
      endif
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
      nva=ncp*nlev*(ntg+1)
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
      mavfound=0
c      write(*,*) 'runlog=',runlog
      open(lun_rlg,file=runlog,status='unknown')
      read(lun_rlg,*,err=888) nlhead,ncol
      do i=2,nlhead
         read(lun_rlg,*)
      end do
888   continue  !  Continue to support old format runlogs
      if(debug) write(*,*)' Main loop...',nspectra,mspectra
      do ispec=1,mspectra         !  Main fitting loop over spectra
141     call read_runlog(lun_rlg,col1,specname,iyr,iset,zpdtim,
     &  oblat,oblon,obalt,asza,zenoff,azim,osds,
     &  opd,fovi,fovo,amal,ifirst,ilast,graw,possp,bytepw,zoff,snr,apf,
     &  tins,pins,hins,tout,pout,hout,
     &  sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)
        if(debug) write(*,*) 'runlab, istat=',specname, istat
        if(istat.ne.0) exit
        if(col1.eq.':') cycle

        lr=lnbc(specname)
        specname=specname(:lr)
        if(lr.gt.lrmax) lrmax=lr  ! Longest spectrum name
c
        read(lun_ray,'(a,200f11.5)') specray,rdum,rdum,sza_ray,bend,
     &  rdum,zmin,(splos(j),j=1,nlev)
c        write(37,*)zmin,zminwas,sza_ray,sza_raywas,
c     &  (zmin-zminwas)/(sza_ray-sza_raywas)
        zminwas=zmin
        sza_raywas=sza_ray
        specray=specray(:lr)
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
c  Apply air-to-vacuum  & FOV corrections
        if(kspflag.lt.2) graw=graw*riair(lasf,tins,pins,hins)/
     &  riair(frqcen,tins,pins,hins)
c        graw=graw*(1.D0+(amal**2+fovi**2)/16)  ! FOV correction

        resn=0.5d0/opd
        if(resn.gt.resmax) then
           write(*,*) 'Warning: increase parameter resmax: ',specname
           resn=resmax
        endif
        if(1.0001*resn.lt.graw) then
           write(*,*) 'resn,graw=',resn,graw
           stop ' resn < graw'
        endif
c  Measured spectrum must be wider than fitting interval to
c  allow convolution with apodizing/interpolating ILS
        dd=nscycle*resn ! half-width of the slit function in cm-1 (always +ve)
        vbar=0.5d0*graw*(ilast+ifirst)
        hwid=0.5d0*dabs(graw*(ilast-ifirst))
c        write(*,*)ilast,ifirst,vbar,hwid,dd
        if(kspflag.eq.2) dd=0.0
        if(debug)write(*,*)kspflag,nus,nue,vbar,hwid,vbar-hwid,vbar+hwid
        if( nus-dd .lt. vbar-hwid ) kspflag=1   ! Lower window limit < disk file
        if( nue+dd .gt. vbar+hwid ) kspflag=1   ! Upper window limit > disk file
c        if(nint((nus-dd)/graw).lt.ifirst .or.
c     &  nint((nue+dd)/graw).gt.ilast) kspflag=1
c========================================================================
c  Read model & vmr information (SUNRUN.MAV)
c      write(*,*)'specname, runlabmav=',specname, runlabmav
      lr=lnbc(specname)
      if(ncall.ge.1) then
      if(specname.eq.runlabmav) then
         call read_mavfile_body(lun_mav,nlev,nspexi,z,t,p,d,vmr)
         read(lun_mav,'(a)',end=66) mavstring
         if(mavstring(1:14).eq.'Next Spectrum:') then
            read(mavstring(15:),'(a)') runlabmav
         else
            write(*,*) mavstring
            write(*,*)'Failed to find Next Spectrum string'
            stop
         endif
66       continue
         if(index(winfo,' sa_temp ').gt.0) call vadd(t,1,5.,0,t,1,nlev) 
         if(index(winfo,' sa_pres ').gt.0) call vmul(p,1,.95,0,p,1,nlev)

c  Pre-compute absorption coefficient
c         nva=ncp*nlev*(ntg+1)
         call vmov(zero,0,vac,1,nva)
c         write(*,*)' Calling abscoj...', kspflag
         call abscoj(specname,nlev,t,p,d,nspexi,targmol,vmr,vpf,
     &   linefiles,parfile,fzero,grid,ncp,vac,vac(nva+1))
c         write(*,*)' Called abscoj...'

         write(oformat,'(a5,i2.2,a3)')'(1x,a',lrmax+1,',a)'
         write(6,oformat)' Spectrum            ',
     &   colabel(:lcl)//'   AM      OVC        VSF   VSF_err'
         if(mavfound.eq.0) then
            do jtg=1,ntg
               if(speci(jtg).eq.0) then
                 write(*,*) jtg,pars(jtg),'unrecognized or duplicated'
                 stop
               else
                 colabel=colabel(:lnbc(colabel))//'  AM_'//pars(jtg)(:6)
     &           //' OVC_'//pars(jtg)(:6)//' VSF_'//pars(jtg)(:6)
     &           //' VSF_'//pars(jtg)(:lnbc(pars(jtg)))//'_error'
               endif
            end do
            write(lun_col,oformat)' Spectrum            ',
     &      colabel(:lnbc(colabel))
         endif  !  mavfound.eq.0
         mavfound=mavfound+1
      else
         if(mavfound .eq. 0) then
            write(*,'(a)') 'Mavfile data-block is missing: '
            write(*,'(a)') 'Current Spectrum:',specname
            write(*,'(a)') 'Current  Mavfile:',runlabmav
            stop
         endif
      endif          ! (specname.eq.runlabmav)
      endif          ! (ncall.ge.1) then
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
c  FIND the spectral file, return the PATH to the spectrum
      call vmov(zero,0,obsrvd,1,mmp)
c      write(*,*)'kspflag=',kspflag
      if(kspflag.eq.0) then
        call jetspe(specpath,opd,graw,ifirst,ilast,possp,bytepw,nus,
     &  nue,apo_m,interp,zero,zero,
     &  obsrvd,mmp,nmp,startm,gint,rc)
c        write(*,*)'obsrvd=',obsrvd(1),obsrvd(2),obsrvd(3)
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

      if(lsnd.gt.0) then
        read(winfo(lsnd+8:),*) fsn ! =5 for NOMAD
        detnoise=detnoise/fsn  !  detnoise=0.0002 ! NOMAD detector noise
        sphnoise=sphnoise/fsn  !  sphnoise=0.0005 ! NOMAD source photon noise
      endif

c  Add detector noise and source photon noise in quadrature
      rnoise=SQRT(tbar*sphnoise**2+detnoise**2)
      write(*,*)'tbar,rnoise,SNR =',tbar,rnoise,tbar/rnoise
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
      corrld=sqrt(rect**2+((1.+0.4*apo_c)*resn)**2)/gint
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
         do ilev=1,nlev
            call vsma(vac(jva),1,ckm2cm*splos(ilev),spxv(jsp),1,
     &      spxv(jsp),1,ncp)
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
         if(reclen_solarll.eq.101) then
           call solar_pseudo_trans_spec(lun_sts,solarll,
     &     fzero*(1.d0+osds*1.E-6),grid*(1.d0+osds*1.E-6),frac,spts,ncp)
         else
           call solar_pts(lun_sts,solarll,
     &     fzero*(1.d0+osds*1.E-6),grid*(1.d0+osds*1.E-6),frac,spts,ncp)
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
      if (nlev.le.ncell+1) then
         spver(1)=splos(1)
         spver(nlev)=splos(nlev)
      else
         call compute_vertical_paths(ncell,zmin,z,d,spver,nlev)
      endif  ! nlev.le.2
c===================================================================
c Compute vertical and LOS slant columns for the target gases
c Do in reverse order so that the concentration profile of the
c first target gas is retained in CP at the end.
      do jtg=ntg,1,-1
         jspeci=speci(jtg)
         call vmul(vmr(jspeci),nspexi,d,1,cp,1,nlev)
         call vdot(cp,1,spver,1,overcol(jtg),nlev)
         call vdot(cp,1,splos,1,oloscol(jtg),nlev)
         overcol(jtg)=overcol(jtg)*ckm2cm
         oloscol(jtg)=oloscol(jtg)*ckm2cm
      end do
c
      if(ntg.eq.0) then
c   compute total column amounts for air (vmr=1.000)
         call vdot(d,1,spver,1,overcol(1),nlev)
         call vdot(d,1,splos,1,oloscol(1),nlev)
      endif
c=====================================================================
c  Read the initial values of the parameters to be fitted and their estimated
c  a priori variances.
      open(lun_apx,file=ap_file,status='old')
      read(lun_apx,*)
      read(lun_apx,*) apx(n1),apu(n1)   ! Continuum Level
      read(lun_apx,*) apx(n2),apu(n2)   ! Continuum Tilt
      read(lun_apx,*) apx(n3),apu(n3)   ! Continuum Curvature
      read(lun_apx,*) apx(n4),apu(n4)   ! Frequency Shift
      read(lun_apx,*) apx(n5),apu(n5)   ! Zero Offset
      read(lun_apx,*)
      if(ntg.ge.1) read(lun_apx,*) apx(1),apu(1)
      if(ntg.ge.2) read(lun_apx,*) apx(2),apu(2)
      close(lun_apx)
      call vmov(apx(2),0,apx(3),1,ntg-2) ! absorber amount
      apx(n5)=apx(n5)+sngl(zoff)
      call vmov(apu(2),0,apu(3),1,ntg-2)  ! vmr of non-target gases
      apu(n4)=apu(n4)*sngl(1.d0+rdec) ! shift   up from 0.5 4-DEC-98 gct
c==========================================================================
c      write(*,*)' ntg=',ntg
c      write(*,*)'sl: apx(n2)=',n2,apx(n2)
      call vmov(apx,1,cx,1,nfp) ! Initialize to A PRIORI each spectrum
      if(debug)
     &  write(6,*)  'It   CL       CT    CC    FS    ZOFF  RMS/CL'//
     &'  Vfact1  Vfact2  Vfact3  Vfact4  Vfact5  Vfact6'

      ifm=1
      if(index(winfo,' prof_ret_only ').gt.0) ifm=4  ! Profile Retrieval only
      if(index(winfo,' prof_ret+scale ').gt.0) ifm=5  ! Profile Retrieval
      if(index(winfo,'nfov=').gt.0) ifm=3
      if(ifm.eq.1) then   ! FOV center ray only
         call do_retrieval(obsrvd,nmp,apx,apu,slit,nii,
     &    ldec,rdec,spts,spxv,vac,splos,nlev,ncp,ntg,nfp,snr,
     &    corrld,sssss,winfo,debug,mit,nit,calcul,rms,cx,ex,pd,ssnmp)
c      elseif(ifm.eq.2) then  !  Fast approximate FOV integration
c         call do_retrieval2(obsrvd,nmp,apx,apu,slit,nii,
c     &  ldec,rdec,spts,spxv,vac,splos,nlev,ncp,ntg,snr,
c     &  corrld,sssss,winfo,debug,mit,nit,calcul,rms,cx,ex,slpd,pd,ssnmp)
      elseif(ifm.eq.3) then  ! Full numerical FOV integration (slow)
         solzen=sngl(asza+zenoff)
         roc=6378.
         fbar=(ifirst+ilast)*graw/2
         call do_retrieval3(obsrvd,nmp,apx,apu,slit,nii,
     &    z,t,p,solzen,fovo,roc,obalt,wavtkr,fbar,
     &    ldec,rdec,spts,spxv,vac,splos,nlev,ncp,ntg,nfp,snr,
     &    corrld,sssss,winfo,debug,mit,nit,calcul,rms,cx,ex,pd,ssnmp)
      endif

c  Set error bars very small for air, otherwise it will
c  dominate the CO2/Air ratio uncertainties.
      if(index(winfo,' air ').gt.0) ex(1)=1.0e-08

c  Write .spt file for the first MSPT spectral fits
c       write(*,*)'sptpath=',sptpath, ispec,mspt,100*rms/cl
      if(ispec .le. mspt) then
         sptpath=sptfile(:lf-1)//specname
         call write_spt(lun_spt,winfo,sptpath,
     &   obsrvd,calcul,cx,ex,startm+gint*(cx(n4)),
     &   osds*1.0E-06*freq_flag,gint,overcol,pars,asza+zenoff,obalt,
     &   zmin,abs(100*rms/cx(n1)),frac,pd,ssnmp,nmp,nmpfp,ntg,nfp)
      endif

c  Solar_Gas_shift is the difference of the solar shift and the gas shift
c  expressed in terms of the observed spectral point spacing
c  Multiply by GINT to convert to cm-1.
c  Divide by frequency to normalize to a dimensionless stretch which
c  should then  be the same for all windows. Multiply by 10^6 for ppm.
      if(index(winfo,' so ').gt.0) then
         sgshift=1E6*2*gint/frqcen*
     &   solar_gas_shift(cx(n1),cx(n2),obsrvd,calcul,ssnmp,nmp)
c      write(*,*)cx(n1),cx(n2),cx(n4),obsrvd(1),calcul(1),pd(1+n1*nmpfp)
      else
         sgshift=0.0
      endif

c  Write the state vector and uncertainties to the .col file
      call write_col(lun_col,colfile_format,specname,lrmax,nit,rms,
     & sgshift,gint,zmin,oloscol,overcol,cx,ex,ntg,nfp)
c  And to the screen
      call write_col(6,colfile_format,specname,lrmax,nit,rms,
     & sgshift,gint,zmin,oloscol,overcol,cx,ex,ntg,nfp)

      if(debug .and. nit.eq.mit+1) stop 'nit=mit+1' 
c  Output PD's (weighting functions), for subsequent use
c  in deriving averaging kernels.
         if( index(winfo,' ak ') .gt. 0) then
           akpath=akfile(:lnbc(akfile))//specname
           open(lun_ak,file=akpath,status='unknown')
           call fm(lun_ak,winfo,slit,nii,ldec,spts,spxv,
     &     vac,splos,nlev,ncp,rdec,sssss,cx,ntg,nfp,calcul,pd,nmp)
           ynoise=2.5*cx(n1)*corrld/sngl(0.1d0+snr)
c  Skip levels representing the cells, so start ilev at ncell+1.
           write(lun_ak,*)(splos(ilev)*cp(ilev),ilev=ncell+1,nlev)
           write(lun_ak,*) pout/1013.25
           write(lun_ak,*)(p(ilev),ilev=ncell+1,nlev)
           write(lun_ak,*)(ynoise/apu(ifp),ifp=1,nfp)
           write(lun_ak,*)((apx(ifp)-cx(ifp))*ynoise/apu(ifp),ifp=1,nfp)
           close(lun_ak)
         endif

      totit=totit+nit
      avgrms=avgrms+1
      avgcl=avgcl+abs(cx(n1)/rms)
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
         write(6,*)' Total number of mavblocks found =',mavfound
         write(6,*)' Average % RMS fit =',100*avgrms/avgcl
         avgtc=tottc/toterr
         rmstc=sqrt(abs(tottc2/toterr-(avgtc-1.)**2))
         write(6,*)' Average Totcon =',avgtc,' +- ',rmstc
      endif
      close(lun_rlg)
      close(lun_ray)
      close(lun_mav)
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
