c  Spectral Fitting Program GFIT (see GGG.DOC for description of latest changes)
      implicit none
c
      integer*4 rc,lunq,lunr,lunt,lunv,luna,lunb,lunm,lcolon,
     & i,j,k,i0,i1,kn2,kn3,nop,ls,ncol,iseed,jj,n1,n2,n3,n4,istat,
     & lnbc,ierr,iyr,iset,nsh,krank,nmpfp,jsp,jpd,ir,nn,lp,lc,lr,
     & kconv,nconv,         ! Indicates status of convergence (2, 1, 0 )
     & kcp1, kcp2,          ! first and last primative grid indices
     & ifirst,              ! index of first spectral point in disk file
     & ilast,               ! index of last spectral point in disk file
     & apo_m,apo_c,defapo,  ! apodization (0,1,2,3,4) employed for spectral fits
     & ldec,                ! factor by which slit function is oversampled wrt GRID
     & nscycle,             ! slit function breadth (in cycles of sin(x)/x )
     & mmode,nmode,         ! Number of vibrational modes
     & kspecflag,           ! =0 if specified interval of spectrum is readable
     & possp,               ! Length of attached header in bytes
     & bytepw,              ! Bytes per data word, usually 2 (I*2) or 4 (R*4)
     & interp,              ! interpolation in measured spectrum (constant)
     & nlhead_ggg,          ! number of header lines in .ggg file
     & freq_flag,          ! =1  presents spectral fits in solar rest frame.
                            ! =0  presents spectral fits in atmosphere rest frame.
     & mavfound,            ! number of .mav data-blocks found
     & mavwithout,          ! number of .mav data-blocks not found
     & mit, nit,            ! number of iterations allowed to achieve convergence
     & totit,               ! total number of iterations for all the spectra
     & mmp, nmp, imp,       ! number of measured spectral points
     & nhw, nhwmax,         ! Half-width (in grid points) of ILS
     & mii, nii, niimax,    ! Dimension of interpolated (by factor LDEC) ILS
     & kgas,kiso,
     & ncp,                 ! number of pre-calculated VACs
     & nspec,ispec,         ! number of spectra that will be fitted
     & nlevray,             ! number of levels in .ray file
     & mlev,nlev,klev,jlev, ! maximum permitted number of levels
     & mtg, ntg, jtg,       ! number of target molecules
     & mfp, nfp, kfp, jfp,  ! Total number of fitted parameters (usually NTG+4)
     & mva,nva,jva,         ! maximum number of absorption coefficients
     & msp,                 ! maximum dimension of SPXV array 
     & mspt,                ! Maximum number of .spt files written
     & mspexi,nspexi,jspexi ! number of different species listed in PARFILE
c
c  Executable size is approximately  4 + 4.5*MVA  Mbytes
c  MVA = 41200000 produces 204 Mbyte executable
c  MVA = 35000000 produces 170 Mbyte executable
c  MVA = 20000000 produces  84 Mbyte executable
      parameter (lunr=14,lunt=16,lunq=17,lunv=18,luna=20,lunb=21,
     & lunm=22,
     & nscycle=25,mspexi=140,mlev=151,mmp=360000,i0=0,i1=1,
     & mva=60000000,msp=mva/8,mii=103847,mtg=10,mfp=mtg+4,
     & mmode=30,mspt=2500)

      integer*4
     & ip(mfp),         ! used by SHFTI & SCOV2
     & targmol(mspexi), ! Group assignment for each specie.
     & speci(mtg)       ! speci # of the parent isotopolog of each target gas
c
      real*4 
     & tau,dd,zero,wzo,
     & unity,rdum,sumr2,big,tiny,thresh,cl,wshift,wtilt,wcl,wcle,
     & wrms,xo,logrp,xl,dz,fs,cntuum,dxlimit,avgrms,avgcl,
     & gasdev,fr,xfr,
     & cfamp,cffreq,cfphase,tot,cc,cr,misalign,
     & tns,tds,tng,tdg,sgshift,wsgshift,gs,gg, residual      

      real*4
     & tpd(mfp),         !
     & ckm2cm,           ! convert km to cm
     & tcbar,            ! The average transmittance over the window.
     & resoln,           ! desired spectral resolution (FETSPE input parameter)
     & corrld,           ! Factor=[apodixed resolution]/[spectral point spacing]
     & ynoise,           ! Measurement Noise = CORRLD*CX(NTG+1)/SNR
     & foff,             ! frequency offset by which spectrum is to be resampled
     & obsrvd(mmp),      ! observed spectrum returned by FETSPE
     & calcul(mmp),      ! calculated spectrum after convolution/decimation
     & curv(mmp),        ! 
     & noise(mmp),       ! 
     & resids(mmp+mfp),  ! residual spectrum after scaling and tilting
     & rms,rwas,         ! RMS of the spectral fit 
     & var,              ! variance = (rms*corrld)**2 * NMP/(NMP-NTG-3)
     & pd((mmp+mfp)*mfp),! partial differentials (after convolution)
     & vac(mva),         ! volume absorption coefficients (VAC)
     & spxv(msp),        ! SPXV(i,j) = SUM [ SP(k).VAC(i,k,j) ]
     & zmin,             ! lowest altitude along central ray path 
     & slit(mii),        ! slit function
     & vmr(mlev*mspexi), ! buffer for vmr's
     & vpf(mlev*mspexi), ! buffer for vibrational partition functions
c     & conc(mlev,mtg),   ! concentration profiles of target gases
     & cp(mlev),         ! concentration profiles of target gas
     & ax(mfp),          ! "a priori" values of retrieved parameters
     & apu(mfp),         ! "a priori" uncertainties of variables
     & cx(mfp),          ! current values of variables
     & dx(mfp),          ! adjustment to values of variables
c     & dxwas(mfp),       ! last adjustment to values of variables
     & ex(mfp),          ! estimated uncertainty in variables
     & wk(mfp),          ! workspace used by SHFTI subroutine and others
     & overcol(mtg),     ! original (unscaled) vertical column abundance
     & oloscol(mtg),     ! original (unscaled) line-of-sight column abundance
     & airmass(mtg),     ! Ratio: oloscol/overcol
     & splos(mlev),      ! array of line-of-sight slant paths
     & spver(mlev),      ! array of vertical slant paths
     & z(mlev),          ! altitudes of levels (km)
     & t(mlev),          ! temperatures of levels (K)
     & p(mlev),          ! pressures of levels (atm.)
     & d(mlev),          ! number densities at levels (molec.cm-3)
     & clval, clerr,     ! initial continuum level and its uncertainty
     & ctval, cterr,     ! initial continuum tilt and its uncertainty
     & fsval, fserr,     ! initial frequency shift and its uncertainty
     & zoval, zoerr,     ! initial zero offset and its uncertainty
c                          Note that zoval affects Vfact not Verr
c                          Note that zoerr affects Verr but not vfact
     & soval, soerr,     ! solar spectrum scale factor
     & tgval, tgerr,     ! initial scale factor for the first target gas
     & ntval, nterr      ! initial scale factor for the next target gases
      parameter (zero=0.,unity=1.,tau=6.e-06,big=1.e+18,tiny=1.e-36,
     & ckm2cm=1.0E+05)
c
      real*8
     & tottc,tottc2,toterr,avgtc,rmstc,riair,
     & frqcen,           ! central frequency (cm-1) of spectral window
     & width,            ! width (cm-1) of spectral window
     & fzero,            ! frequency of the zero'th computed VAC
     & startm,           ! frequency of first returned measured point
     & nus,              ! microwindow starting frequency (cm-1)
     & nue,              ! microwindow ending frequency (cm-1)
     & graw,             ! spacing of raw spectrum (cm-1) from GETINFO
     & gint,             ! spacing of returned spectrum = graw/interp
     & grid              ! spacing of primitive spectrum 

      real*8 
     & oblat,           ! observation latitude (deg).
     & oblon,           ! observation longitude (deg).
     & obalt,             ! observation altitude (km)
     & zpdtim,           ! Time of ZPD (UT hours)
     & asza,             ! astronomical solar zenith angle (unrefracted)
     & zenoff,           ! zenith pointing offset
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
     & sia,              ! Solar Intensity (Average)
     & sis,              ! Solar Intensity (SD)
     & aipl,             ! Airmass-Independent Path Length (km)
     & lasf,             ! laser frequency (e.g. 15798.03 cm-1)
     & wavtkr,           ! suntracker operating frequency (e.g. 9900 cm-1)
     & eorv,             ! Earth-Object Radial velocity (m/s)
     & ervc,             ! Earth Rotational Velocity Component (m/s)
     & opd,              ! Optical path difference (cm) of interferogram
     & resn,             ! 0.5d0/opd = half width of SINC function in cm-1
     & resmax,           ! maximum value of resn
     & rect,             ! frqcen*(fovi**2+amal**2)/8 = width of rectangle(cm-1)
     & frac,             ! fractional size of FOVO compared with solar diameter
     & dopp,             ! Earth-SUN Doppler stretch.
     & resnog,           ! RESN / GRID
     & rectog,           ! RECT / GRID
     & rdec,             ! Ratio: GINT/GRID = spectral/primitive point spacings
     & shift,            ! frequency shift (cm-1)
     & sh,               ! frequency shift (cm-1)
     & ddum              ! double-precision dummy variable
c
      character
     & colfile*80,       ! output file containing column amounts
     & akfile*80,        ! output file of averaging kernels
     & sptfile*80,       ! output file of ascii spectral fits
     & mavfile*80,       ! file containing T/P & VMR at user-chosen levels
     & rayfile*80,       ! file of slant paths at user-chosen levels
     & parfile*80,       ! path to "molparam.dat"
     & apvalerr*80,      ! path to a priori variable values and uncertainties
     & llsize*80,        ! path to llsize.dat file
     & linefiles*400,    ! paths to linelists
     & solarll*80,       ! solar linelist
     & version*64,       ! version number
     & gsversion*64,     ! GSETUP version number
     & dplist*80,        ! Data Partition List (e.g. m4part.lst)
     & runlog*80,        ! name of occultation file
     & runlab*21,        ! name of spectrum read from runlog file
     & runlabmav*21,     ! name of spectrum read from mavfile
     & mavstring*48,
     & string*48,
     & apf*2,            ! apodization function (e.g. BX, TR)
     & gasname*8,! names of species in PARFILE
c     & tname*8,  ! names of species in PARFILE prefixed by "t"
c     & fullname*10,! full names of species in PARFILE
     & winfo*128,         !  window information (command line)
     & pars(mtg)*9,     ! parameters to fit
     & specray*21,       ! name of spectrum read from .ray file
     & specpath*128,     ! location of spectrum
c     & fdate*26,         ! current date
c     & getlog*8,         ! investigator
c     & hostnam*12,       ! computer hostname
     & levfile*80,       ! name of file containing fitting levels
     & amodel*80,        ! name of atmospheric model
     & vmrset*80,        ! name of vmr set
     & window*80,        ! name of window list
     & col1*1,           ! first column of runlog record
     & colabel*500       ! column labels for .col file
c
      logical debug,cf
c
      integer*4
     &   dgen(mmode),  ! degeneracy of vibrational modes
     &   molewt     ! Molar Mass
c
      real*4
     &   tdrpf,        ! T-Dependence of Rotational Partition Function
     &   fia,delta,epsilon,  ! Fractional Isotopic Abundance
     &   vibfrq(mmode) ! Array of vibrational frequencies
c
      character
     &   speci_id*24
c
      data speci/mtg*0/
      data version/' GFIT      Version 2.40.2    5-Dec-2006    GCT  '/

      write(6,*)
      write(6,88)version
      nconv=2
      iseed=44444
c
c  Read runlog, model, vmrset & window information from input file (.ggg)

c      go to 10
      read(5,*)nlhead_ggg
      read(5,88)gsversion
      read(5,88)dplist
      read(5,88)apvalerr
      read(5,88)runlog
      read(5,88)levfile
      read(5,88)amodel
      read(5,88)vmrset
      read(5,88)mavfile
      read(5,88)rayfile
      read(5,88)parfile
      read(5,88)window
      read(5,88)llsize
      read(5,88)linefiles
      read(5,88)solarll
      read(5,88)akfile
      read(5,88)sptfile
      read(5,88)colfile
      read(5,88)winfo
      go to 20

c     This section only for running debugger in Digital Fortran
 10   open(62,file="n2o_2439.uow00224.ggg",status="old")
      read(62,*)nlhead_ggg
      read(62,88)gsversion
      read(62,88)dplist
      read(62,88)apvalerr
      read(62,88)runlog
      read(62,88)levfile
      read(62,88)amodel
      read(62,88)vmrset
      read(62,88)mavfile
      read(62,88)rayfile
      read(62,88)parfile
      read(62,88)window
      read(62,88)llsize
      read(62,88)linefiles
      read(62,88)solarll
      read(62,88)akfile
      read(62,88)sptfile
      read(62,88)colfile
      read(62,88)winfo
      close (62)

 20   continue
      lc=index(winfo,'#')
      if(lc.gt.0) winfo=winfo(:lc-1)
      read(winfo,*) frqcen,width,mit,defapo,interp,freq_flag
      if (interp.le.0) stop 'interp must be > 0'
      lp=lnbc(winfo)
      lc=index(winfo,':')
      write(6,88)winfo(:lp)
      call lowercase(winfo)
      call substr(winfo(lc+1:),pars,mtg,ntg)
      if(ntg.gt.mtg) write(*,*)' Warning: Increase MTG to',ntg
      if(debug) write(*,*)'ntg1=',ntg
c
      if( index(winfo,'debug') .gt. 0 ) then
         debug=.true.
      else
         debug=.false.
      endif
c
      cf=.false.
      if( index(winfo,' cf ') .gt. 0 ) then
         nconv=3
         write(*,*)'Fitting channel fringes'
         cf=.true.
      endif
      nus=frqcen-width/2
      nue=frqcen+width/2
      if(debug) write(*,*)'nus,nue=',nus,nue
c
c  Read in names of isotopomers
      if(debug) write(*,*) parfile
      open(lunb,file=parfile,status='old')
      do jspexi=1,mspexi
         call read_isotop(lunb,kgas,kiso,gasname,speci_id,
     &   fia,delta,epsilon,molewt,tdrpf,vibfrq,dgen,nmode,mmode,istat)
         if(istat.ne.0) go to 77
         call lowercase(gasname)
         targmol(jspexi)=0
         do jtg=1,ntg
            if(char(kiso+48)//gasname.eq.pars(jtg)) targmol(jspexi)=jtg
            if( gasname//' '.eq.pars(jtg) .or. 
     &          't'//gasname.eq.pars(jtg)) then
                if(targmol(jspexi).eq.0) targmol(jspexi)=jtg
            endif
         end do
      end do
      read(lunb,*,end=77)
      stop ' The number of species listed in PARFILE exceeds MSPECI'
77    close(lunb)
      nspexi=jspexi-1
      if(debug) write(*,*)'targmol=',targmol
c
      do jspexi=nspexi,1,-1
           if(targmol(jspexi).gt.0) speci(targmol(jspexi))=jspexi
      end do
      if(debug) write(*,*)'ntg=',ntg
      if(debug) write(*,*)'speci=',speci
c
      colabel=
     &'   Spectrum           Nit  CL  Tilt  Sh  S-G  Zoff  RMS/CL  Zmin'
      do jtg=1,ntg
         if(speci(jtg).eq.0) then
            write(*,*) jtg,pars(jtg),'unrecognized or duplicated'
            stop
         else
            colabel=colabel(:lnbc(colabel))//'   AM_'//pars(jtg)(:4)
     &      //' OVC_'//pars(jtg)(:6)//' VF_'//pars(jtg)(:6)
     &      //' VF_'//pars(jtg)(:lnbc(pars(jtg)))//'_error'
         endif
      end do
c
      nfp=ntg+4
      n1=ntg+1
      n2=ntg+2
      n3=ntg+3
      n4=ntg+4
c=====================================================================
c  Read the initial values of the parameters to be fitted and their estimated
c  a priori variances.
      open(19,file=apvalerr,status='old')
      read(19,*)
      read(19,*) clval,clerr
      read(19,*) ctval,cterr
      read(19,*) fsval,fserr
      read(19,*) zoval,zoerr
      read(19,*) soval,soerr
      read(19,*) tgval,tgerr
      read(19,*) ntval,nterr
      close(19)
c========================================================================
c  Read the spectral headers, compute the widths of their slit functions, 
c  This ensures that the spectral headers are all OK before investing alot
c  of time computing absorption coefficients and starting the retrievals.
 88   format(a)
      resmax=0.0d0
      niimax=0
      avgrms=0.0
      tottc=0.0
      tottc2=0.0
      toterr=0.0
      avgcl=0.0
c      nspec=0
      grid=0.666666d-06*frqcen
c  Get the header/runlog information pertaining to RUNLAB
      if(debug) write(*,*)' Opening .ray file: '//
     & rayfile(:lnbc(rayfile))
      open(lunr,file=rayfile,status='old')
      read(lunr,*)nn,ncol
      nlevray=ncol-7  !  First 7 columns are (spec,Zobs,Pobs,SZA,Bend,FOV,Zmin)
      nlev=nlevray
      if(nlev.gt.mlev) then
         write(*,*) 'Warning: NLEV > MLEV :',nlev,mlev
         stop
      endif
      call skiprec(lunr,nn-1)
      ir=1+index(runlog,'.')
      if(debug) write(*,*)' Opening runlog: '//runlog(:lnbc(runlog))
      open(lunq,file=runlog,status='old')
      read(lunq,*)
      open(lunm,file=mavfile,status='old')
      read(lunm,*)
      read(lunm,'(14x,a)')string
      lcolon=index(string,':')
      read(string(lcolon+1:),'(a)')runlabmav

c  Pre-screening Loop over spectra
c  Checks that the spectra, .mav, and .ray files are all okay,
c before investing alot of time on computing absorption coeffs.
      do ispec=1,9999999
131   call read_runlog(lunq,col1,runlab,iyr,iset,zpdtim,
     & oblat,oblon,obalt,asza,zenoff,opd,fovi,fovo,amal,ifirst,
     & ilast,graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     & tout,pout,hout,lasf,wavtkr,sia,sis,aipl,istat)
c      write(*,*)'called read_runlog, istat,graw=',istat,graw
      if(istat.ne.0) go to 89

      lr=lnbc(runlab)
      runlab=runlab(:lr)
c
c  Check that .ray file is readable and consistent with runlog
c  Note that in gfit, the cell_length is read as the first
c  element of splos, whereas in gsetup it is a separate variable.
      read(lunr,*)specray,rdum,rdum,rdum,rdum,rdum,zmin,
     &(splos(j),j=1,nlevray)
c      write(*,*)nlevray,splos
      specray=specray(:lr)
      if(specray.ne.runlab) then
         write(*,*) specray//'   '//runlab
         stop 'spectrum mismatch 1'
      endif

c  If MIT > 0, check that requested spectrum is on disk,
c  and that it covers the specified spectral interval.
      kspecflag=0
      call gindfile(dplist,runlab,specpath)
      if(lnbc(specpath).eq.0) then
         write(*,*) runlab,' not found on disk'
         kspecflag=2
      endif
      resn=0.5d0/opd
c      if(resn.lt.graw) resn=graw
c  Measured spectrum must be wider than fitting interval to
c  allow convolution with apodizing/interpolating ILS
      dd=nscycle*resn ! half-width of the slit function in cm-1
      dd=0.0d0
c      write(*,*)opd,nus,nue,dd,ifirst,ilast,graw
      if(debug) write(*,'(a21,a13,i5,a4,i5,a5)') runlab,' encompasses ',
     & nint(ifirst*graw),' to ',nint(ilast*graw),' cm-1'
      if(nint((nus-dd)/graw).lt.ifirst .or.
     & nint((nue+dd)/graw).gt.ilast) then
c         write(*,*)' Requested interval extends from ',nus,' to ',nue
c         write(*,*)' Required interval extends from ',nus-dd,
c     & ' to ',nue+dd
c         write(*,*)' Spectrum covers ',ifirst*graw, ilast*graw
         kspecflag=1
      endif
c      write(*,*)'kspecflag= ',kspecflag
      if(mit.gt.0 .and. kspecflag.gt.0) goto 131 ! skip missing/partial spectrum
      if(mit.eq.0 .and. kspecflag.eq.1) goto 131 ! skip missing/partial spectrum
c
      if(dabs(snr).lt.1.d0) write(6,*) 
     & 'Warning: S/N ratio of run '//specray//' =',snr
      if(dabs(opd).lt.1.d0) write(6,*) 
     & 'Warning: OPD of run '//specray//' =',opd
c
      if(resn.gt.resmax)resmax=resn
c
      ldec=1+int(16*grid/resn)
      nhw=nint(nscycle*resn/grid)
      nii=1+2*ldec*nhw
      if(nii.gt.niimax) niimax=nii
c
c      nspec=nspec+1
      end do ! jspe=1,999999  ! Loop over spectra
c-------------------------------------------------------
 89   nspec=ispec-1
      close(lunr)
      close(lunq)
      if(nspec.le.0)
     & write(*,*) 'GFIT Warning: no spectra found for this interval'
c
      if(niimax.gt.mii) then
        write(6,*) 'STOP: increase MII to ',niimax
        stop
      endif
c======================================================================
c  Define spectral limits of VAC calculation, check that arrays are large
c  enough, choose linelists, and finally pre-compute VACs.
12    kcp1=int(nus/grid)
      kcp2=int(nue/grid)
      nsh=int(2+2*resmax/grid)
      nhwmax=nint(nscycle*resmax/grid)
c      write(6,*)'gfit: nus, nue, resmax, nscycle, grid ',nus,nue,resmax,
c     &grid,nscycle
c      write(6,*)'gfit: kcp1, kcp2, nsh, nhwmax ',kcp1,kcp2,nsh,nhwmax
      ncp=kcp2-kcp1+2*nhwmax+2*nsh
      fzero=grid*(kcp1-nsh-nhwmax)
c      write(6,*)'gfit: ncp, fzero ',ncp,fzero
      nva=ncp*nlev*(ntg+1)
      if(debug) write(6,*)'mva, nva =',mva,nva
      write(*,*)' Grid=',grid,'cm-1' 
      write(*,*)' Using ',float(nva+ncp)/mva,' of allocated memory'
      if(nva+ncp.gt.mva) then
        write(6,*)'Warning: loss of accuracy may result from increased'
        write(6,*)'primitive grid spacing forced by insufficient MVA.'
        write(6,*)'To avoid this, Increase MVA from',mva,' to',nva+ncp
        grid=1.0022d0*grid*(nva+ncp)/mva
        go to 12
      endif
c=========================================================================
      if(nfp*ncp.gt.msp) then
         write(6,*)'Increase dimension of MSP from ',msp,' to ',ncp*nfp
         stop
      endif
c
c      call getenv('HOSTNAME',hostnam)
c  Read one measured spectrum at a time, calculate matching synthetic,
c  perform retrieval, and write results to disk (OUTFILE).
      open(lunt,file=colfile,status='unknown')
      write(lunt,*)nlhead_ggg+2,7+4*ntg
      write(lunt,88)version
c     &//fdate()//getlog()//hostnam(:lnbc(hostnam))  ! cpp
      write(lunt,88)gsversion(:lnbc(gsversion))
      write(lunt,88)dplist(:lnbc(dplist))
      write(lunt,88)apvalerr(:lnbc(apvalerr))
      write(lunt,88)runlog(:lnbc(runlog))
      write(lunt,88)levfile(:lnbc(levfile))
      write(lunt,88)amodel(:lnbc(amodel))
      write(lunt,88)vmrset(:lnbc(vmrset))
      write(lunt,88)mavfile
      write(lunt,88)rayfile
      write(lunt,88)parfile
      write(lunt,88)window(:lnbc(window))
      write(lunt,88)llsize(:lnbc(llsize))
      write(lunt,88)linefiles(:lnbc(linefiles))
      write(lunt,88)solarll(:lnbc(solarll))
      write(lunt,88)akfile(:lnbc(akfile))
      write(lunt,88)sptfile(:lnbc(sptfile))
      write(lunt,88)colfile(:lnbc(colfile))
      write(lunt,88)winfo(:lp)
      write(lunt,88) colabel(:lnbc(colabel))
c
      totit=0
      mavfound=0
      mavwithout=0
      open(lunr,file=rayfile,status='old')
      read(lunr,*)nn,ncol
      nlevray=ncol-7  ! First 6 columns are (spec,Zobs,Pobs,SZA,Bend,FOV,Zmin)
      call skiprec(lunr,nn-1)
      open(lunq,file=runlog,status='unknown')
      read(lunq,*)
      if(debug) write(*,*)'nspec=',nspec
      do ispec=1,nspec         !  Main fitting loop over spectra
141     call read_runlog(lunq,col1,runlab,iyr,iset,zpdtim,
     &  oblat,oblon,obalt,asza,zenoff,opd,fovi,fovo,amal,ifirst,
     &  ilast,graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     &  tout,pout,hout,lasf,wavtkr,sia,sis,aipl,istat)
c        write(*,*)'ispec,nspec,istat,runlab=',ispec,nspec,istat,runlab
        if(istat.ne.0) then
            stop 'error in read_runlog'
        endif

c  Apply air-to-vacuum  & FOV corrections
        graw=graw*riair(lasf,tins,pins,hins)/
     &  riair(frqcen,tins,pins,hins)

c        graw=graw*(1.D0+(amal**2+fovi**2)/16)  ! FOV correction

        lr=lnbc(runlab)
        runlab=runlab(:lr)
c
c  Select apodisation function
        if(apf.eq.'BX') then
           apo_m=defapo
           apo_c=defapo
        else  ! if the measured spectra are already apodized
           apo_m=-1  ! for perfect representation of synthetic spectra
           if(apf.eq.'N1') then
             apo_c=1
           elseif(apf.eq.'N2') then
             apo_c=2
           elseif(apf.eq.'N3') then
             apo_c=3
           elseif(apf.eq.'TR') then
             apo_c=4
           else 
             write(6,*)runlab//apf//' ???  Unknown apodization function'
             stop
           endif
        endif
c        apo_m=-1  ! for perfect representation of synthetic spectra
c
        read(lunr,*)specray,rdum,rdum,rdum,rdum,rdum,zmin,
     &  (splos(j),j=1,nlev)
        specray=specray(:lr)
c 34     format(a21,6f10.4,151f8.3)
        if(specray.ne.runlab) then
          write(6,*) specray,runlab
          stop 'spectrum mismatch 2'
        endif
c
c  If MIT > 0, check that requested spectrum is on disk,
c  and that it covers the spexified spectral interval.
        kspecflag=0
        call gindfile(dplist,runlab,specpath)
        if(lnbc(specpath).eq.0) kspecflag=2
        resn=0.5d0/opd
c        if(resn.lt.graw) resn=graw
        dd=nscycle*resn ! half-width of the slit function in cm-1
        dd=0.0
        if(nint((nus-dd)/graw).lt.ifirst .or.
     &   nint((nue+dd)/graw).gt.ilast) kspecflag=1
        ls=lnbc(runlog)
c========================================================================
c  Read model & vmr information (SUNRUN.MAV)
c      write(*,*)'runlab, runlabmav=',runlab, runlabmav
      if(runlab.eq.runlabmav) then
           mavfound=mavfound+1
           call read_mav(lunm,mlev,nlev,nspexi,z,t,p,d,vmr)
           if(debug) write(*,*)'z,t,p,d=',z(1),t(1),p(1),d(1)
           if(debug) write(*,*)'z,t,p,d=',z(2),t(2),p(2),d(2)
           read(lunm,'(a)',end=66) mavstring
           if(mavstring(1:14).eq.'Next Spectrum:') then
              read(mavstring(15:),'(a)') runlabmav
           else
              write(*,*) mavstring
              write(*,*)'Failed to find Next Spectrum string'
              stop
           endif
66         continue
           if(index(winfo,' sa_temp ').gt.0)
     &     call vadd(t,1,5.0,0,t,1,nlev)
           if(index(winfo,' sa_pres ').gt.0)
     &     call vmul(p,1,0.95,0,p,1,nlev)
c
c  Pre-compute absorption coefficient
           if(debug) write(*,*) 'Calling abscoh',fzero,fzero+grid*ncp
           call vmov(zero,i0,vac,i1,nva)
           call abscoh(nlev,t,p,d,nspexi,targmol,vmr,vpf,llsize,
     &     linefiles,parfile,fzero,grid,ncp,vac(1),vac(nva+1))
        call vdot(vac,1,vac,1,sumr2,ncp)
           write(6,88)runlab//' '//
     &     'Nit  CL  Tilt  Sh   S-G  Zoff RMS/CL'//
     &     '  Zmin Airmass Org.VCol      VF    VF_err'
      else
         mavwithout=mavwithout+1
         if(mavfound .eq. 0) then
            write(*,*) 'Mavfile data-block is missing: '
            write(*,*) 'Current Spectrum:',runlab
            write(*,*) 'Current Mavfile:',runlabmav
            stop
         endif
      endif          ! (runlab.eq.runlabmav)
      if(mit.gt.0 .and. kspecflag.gt.0) goto 141 ! skip missing/partial spectrum
      if(mit.eq.0 .and. kspecflag.eq.1) goto 141 ! skip missing/partial spectrum
c=========================================================
      resn=0.5d0/opd
c      if(resn.lt.graw) resn=graw
      if(index(winfo,' sa_fovi ').gt.0) fovi=fovi*1.07
      rect=frqcen*(fovi**2+amal**2)/8  ! old code
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
c  Read spectrum RUNLAB into array OBSRVD over the interval NUS to NUE cm-1
      resoln=0.0    ! use the full resolution of the actual spectrum
      foff=0.0    !  use original measured spectral points (do not resample)
c-------------------------------------------------------------------------
C  FIND the spectral file, return the PATH to the spectrum
      call vmov(zero,0,obsrvd,1,mmp)
      if(kspecflag.eq.0) then
        if(debug) write(*,*)'Calling JETSPE...',apo_m,apo_c
        call jetspe(specpath,opd,graw,ifirst,ilast,possp,bytepw,nus,
     &  nue,apo_m,interp,foff,resoln,slit,mii,nscycle,
     &  obsrvd,mmp,nmp,startm,gint,rc)
        if(rc.ne.0) then
           write(6,*)' Error in IETSPE. Spectrum ',runlab,rc,nmp
           write(6,*)' This error should never happen'
           kspecflag=1
        endif
      endif
      if(kspecflag.eq.2) then
        startm=nus
        gint=graw/interp
        nmp=1+(nue-nus)/gint
      endif
       if(nmp.gt.mmp) stop 'increase parameter MMP'
      if(index(winfo,' sa_snr ').gt.0) then
      do i=1,nmp
      noise(i)=0.0001*gasdev(iseed)*opd ! add 0.25% rms random noise
c      noise(i)=noise(i)*(1+0.01*4*i*(nmp-i)/nmp/nmp)  ! add +- 1% Cont Curv 
c      noise(i)=noise(i)*(1+0.005*obsrvd(i))  !  0.5%  detector non-linearity
      end do
      endif
      iseed=iseed+1113
      nmpfp=nmp+nfp
c==========================================================================
c  Pre-compute ILS (oversampled by a factor LDEC) and normalize each
c  interleaved component to unity.
      rdec=gint/grid
      nhw=nint(nscycle*resn/grid) 
      ldec=1+int(16*grid/resn)
      nii=1+2*nhw*ldec
      resnog=ldec*resn/grid
      rectog=ldec*rect/grid
      if(debug) write(*,*)'Calling profzl...',apo_c,nii,resnog,rectog
      call profzl(apo_c,nii,resnog,rectog,0.0d0,slit)
c      write(*,*)resn, rect, ldec, grid
c      write(*,*)resnog, rectog, nhw, nii
      do j=ldec,1,-1
        call vdot(slit(j),ldec,unity,0,tot,2*nhw)
        call vmul(slit(j),ldec,1.0/tot,0,slit(j),ldec,2*nhw)
      end do
      slit(nii)=slit(nii)/tot
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
         do jlev=1,nlev
            call vsma(vac(jva),i1,ckm2cm*splos(jlev),spxv(jsp),i1,
     &      spxv(jsp),i1,ncp)
c            call vdot(vac(jva),1,vac(jva),1,sumr2,ncp)
c            if(jlev.eq.1 .and. jtg.eq.1)
c     &      write(*,*)'GFIT: jtg,jlev,jva,sumr2=',jtg,jlev,jva,sumr2
            jva=jva+ncp
         end do
         jsp=jsp+ncp
      end do
c
      dopp=0.0d0
      kn2=1+ncp*(ntg+1)   ! Address of start of solar spectrum  
c   Add solar optical thickness spectrum to non-target ones in SPXV(JSP)
c   Use SPXV(JSP+NCP) as work space for Voigt functions.
      if( index(winfo,' so ') .gt. 0) then
         if    (runlab .eq. 'phg92258.600') then
            dopp=-0.04d-5
         elseif(runlab .eq. 'pin92258.600') then
            dopp=-0.22d-5
         elseif(runlab .eq. 'psl3.sun.1') then
            dopp=+2.68d-5
         elseif(runlab .eq. 'psl3.sun.2') then
            dopp=+2.42d-5
         elseif(runlab .eq. 'psl3.sun.3') then
            dopp=+2.37d-5
         elseif(runlab .eq. 'psl3.sun.4') then
            dopp=+2.3d-5
         elseif(runlab .eq. 'pat3.f09ss.sun') then
            dopp=-0.34d-5
         elseif(runlab .eq. 'pat3.f12ss.sun') then
            dopp=+0.23d-5
         elseif(runlab .eq. 'pat3.f03ss.sun') then
            dopp=-0.13d-5
         elseif(runlab .eq. 'pat3.f04ss.sun') then
            dopp=+0.89d-5
         elseif(runlab .eq. 'camy-peyret.solar_ir') then
            dopp=-0.06d-5
         elseif(runlab .eq. '901218R0.003') then
            dopp=-1.11d-6
         elseif(runlab .eq. '901218R0.006') then
            dopp=+0.72d-6
         elseif(runlab .eq. '810509R0.004') then
            dopp=+1.55d-6
         elseif(runlab .eq. '810509R0.005') then
            dopp=+2.75d-6
         elseif(runlab(:9) .eq. 'phg06000.')then
              dopp=0.0d0
         elseif(runlab(:9) .eq. 'phg07000.')then
              dopp=0.0d0
         elseif(runlab(:9) .eq. 'phg06000.')then
              dopp=0.0d0
         elseif(runlab(:9) .eq. 'phg07000.')then
              dopp=0.0d0
         else
c compute Earth - Sun doppler stretch
            if(debug) write(*,*)'Calling ZENAZ'
            call zenaz(2,oblat,oblon,obalt,iyr,1,iset,
     &      zpdtim/24.0d0,ddum,ddum,eorv,ervc,
     &      ddum,ddum,ddum)
            dopp=(eorv+ervc)/3.d+08
         endif

         frac=fovo/9.2e-3  ! the sun is 9.2 mrad in diameter on average
         if(debug) write(*,*)'Calling SOLAR_SPECTRUM',grid,frac,dopp
         call solar_spectrum(solarll,fzero*(1.d0+dopp),
     &   grid*(1.d0+dopp),frac,spxv(kn2),ncp)
         if(debug) write(*,*)'Called SOLAR_SPECTRUM'
      else
         call vmov(unity,0,spxv(kn2),1,ncp)
      endif
      if(index(winfo,' sa_zoff ').gt.0) zoff=zoff+0.002
c=========================================================================
c  Error estimation naively assumes spectral points are linearly independent.
c  CORRLD is an estimate of how dependent neighbouring points really are.
c  CORRLD = 1 means that points are truly linearly independent.
c  Variances must be scaled by CORRLD to ensure that when the same spectrum is
c  analyzed under different values of APO or INTERP, the same size errors result
      corrld=sqrt(rect**2+((1.+0.4*apo_c)*resn)**2)/gint
c  Set inverse uncertainties (RAE) to "a priori" values
      apu(1)=tgerr     ! target   increased from 0  6-19-95  GCT
      call vmov(nterr,i0,apu(2),i1,ntg-1)  ! vmr of non-target gases
      apu(n1)=clerr ! cont     increased from 0 7-2-94 GCT
      apu(n2)=cterr ! tilt     increased from 0.1 on 29-Dec-97 GCT
      apu(n3)=fserr*sngl(1.d0+rdec) ! shift   up from 0.5 4-DEC-98 gct
      apu(n4)=zoerr     ! zero offset error
c  Set AX to A PRIORI estimates 
      ax(1)=tgval
      call vmov(ntval,i0,ax(2),i1,ntg-1) ! absorber amount
      ax(n1)= clval             ! continuum level
      ax(n2)= ctval             ! continuum tilt
      ax(n3)=fsval+sngl((startm/grid-kcp1+nsh-1+nhwmax-nhw)/rdec)
      ax(n4)=zoval+sngl(zoff)

c
      call vmov(ax,i1,cx,i1,nfp) ! Initialize to A PRIORI each spectrum
      if(index(winfo,' cl ').gt.0) then
           cx(n1)=(obsrvd(1)+obsrvd(nmp))/2
           if(cx(n1).lt.1.0) cx(n1)=1.0
      else
           apu(n1)=1.
      endif
      if(debug) write(6,*)  'It  Cntuum    Tilt  Shift  ZOFF  RMS/CL'//
     &'  Vfact1  Vfact2  Vfact3  Vfact4  Vfact5  Vfact6'
c
      shift=gint*dble(cx(n3)-ax(n3))
      rwas=big
      kconv=nconv
      do 65 nit=0,mit              ! Spectral fitting iteration loop
c  Calculate from where (RVA) to start using primitive spectrum (PTSPEC)
c  and prevent our usage from ever exceeding array bounds.
        if(dabs(shift) .gt. 1*(resn+grid)) then
          if(debug) write(6,*)' Limiting Frequency shift'
          shift=dsign(1*(resn+dble(grid)),shift)
          cx(n3)=ax(n3)+sngl(shift/gint)
        endif
c  Limit TILT if it exceeds 1.0
c        if( abs(cx(n2)) .ge. 1.0 ) then
c          if(debug) write(6,*)' Limiting TILT'
c          cx(n2)=sign(1.0,cx(n2))
c        endif
c  Calculate spectrum & PD's
        if(debug) write(*,*)'Calling FM....',cx(1),cx(2),cx(n3)
        call fm(winfo,slit,nii,ldec,spxv,ncp,rdec,cx,ntg,calcul,pd,nmp)
        call vdot(calcul,1,calcul,1,sumr2,nmp)
        if(debug) write(*,*)'rms_calcul=',sqrt(sumr2/nmp)
        if(index(winfo,' sa_snr ').gt.0)
     &  call vadd(calcul,1,noise,1,calcul,1,nmp)
c
c  Calculate spectral curvature (to compute misalignment)
        call vsma(calcul(2),i1,-0.5,calcul,1,curv,1,nmp-1)
        call vsma(calcul(1),i1,-0.5,curv(2),1,curv(2),1,nmp-1)
        curv(1)=0.0
        curv(nmp)=0.0
c==========================================================================
c  Calculate residuals 
        if(kconv.eq.nconv) then ! use logarithmic residuals for 1'st convergence
           do imp=1,nmp
             if(calcul(imp).le.tiny*obsrvd(imp).or.
     &       obsrvd(imp).le.tiny*calcul(imp)) then
                resids(imp)=obsrvd(imp)-calcul(imp)
             else
                resids(imp)=calcul(imp)*log(obsrvd(imp)/calcul(imp))
             endif
           end do   ! imp=1,nmp
        else  !                 use linear residuals for subsequent convergences
           call vsub(obsrvd,i1,calcul,i1,resids,i1,nmp) ! residuals
        endif   !  kconv.eq.nconv
c
c=========================================================================
c Calculate RMS fit
        call vdot(resids,i1,resids,i1,sumr2,nmp)
        rms=sqrt(sumr2/nmp)
        if(debug) write(*,*)'%rms=',100*rms
        if(abs(rms).gt.1.E+15) rms=sign(1.E+15,rms) ! prevents rms=Inf

        if(debug) write(*,*)'nit,ntg,cx=',nit,ntg,cx(n1),cx(n2),
     & 1000*shift,100*cx(n4),100*rms/cx(n1),(cx(jtg),jtg=1,min0(6,ntg))
 82     format(a7,i2,1pe11.4,0pf7.3,f6.3,f7.3,f8.4,9f8.4)
c        if(mit.eq.1) go to 96  ! don't compute PD's on last iteration
        if (rms.gt.9*rwas) then  ! Fit got much worse,
           if(debug) write(6,*)'Retracing step',rms,rwas
           call vmul(dx,1,-0.9,0,dx,1,nfp)
        else
           if(debug) write(6,*)'Continuing',rms,rwas
           thresh=(64*kconv-63)*(rms+0.01*abs(cx(n1)))/100000
           if(abs(rms-rwas).lt.thresh .or. nit+2*kconv-1.gt.mit) then
              kconv=kconv-1
              rwas=big
              if(kconv.ge.1) then
                if(cf) then ! Subtract fitted cosine from OBSRVD
                call vdot(calcul,1,unity,0,tcbar,nmp)
                tcbar=tiny+tcbar/nmp
                call fringes(cfamp,cffreq,cfphase,resids,nmp,mmp)
                if(debug) write(6,'(a,3f11.4)')'Channel Fringes:',
     &          cfamp/cx(n1),2*3.14159*gint*(nmp-1)/cffreq,cfphase
                write(6,'(a,3f11.5)')'Channel Fringes:',
     &          cfamp/cx(n1),2*3.14159*gint*(nmp-1)/cffreq,cfphase
                call vramp(resids,1,nmp)
                call vsma(resids,1,cffreq,cfphase,0,resids,1,nmp)
                call vcos(resids,1,resids,1,nmp)
                call vmul(resids,1,calcul,1,resids,1,nmp)
                call vsma(resids,1,-(cfamp/tcbar),obsrvd,1,obsrvd,1,nmp)
                endif
                call vsub(obsrvd,i1,calcul,i1,resids,i1,nmp) ! residuals
c               nit=nit-1
              endif
           endif ! abs(rms-rwas).lt.thresh .or. nit+2*kconv-1.gt.mit
c
           ynoise=2.5*cx(n1)*corrld/sngl(0.1d0+snr)
c   The last few (i > NMP) elements of RESID & PD contains A PRIORI information
c   RESIDS(nmp+i) holds the values of (AX(i)-CX(i))*RAE(i)*YNOISE
c   while  PD(nmp+i,i) holds RAE(i)*YNOISE, the other elements being zero.
           call vdiv(ynoise,i0,apu,i1,pd(nmp+1),nmpfp+1,nfp)
           call vsub(ax,i1,cx,i1,wk,i1,nfp)
           call vmul(wk,i1,pd(nmp+1),nmpfp+1,resids(nmp+1),i1,nfp)
           if(debug) then
              write(*,122)'ax=',(ax(j),j=n1,n4),(ax(j),j=1,ntg)
              write(*,122)'cx=',(cx(j),j=n1,n4),(cx(j),j=1,ntg)
122           format(a3,7f11.5)
c              do jfp=nmp+1,nmpfp
c              write(6,'(7e11.4)')(pd(jfp+j*nmpfp),j=0,nfp-1),resids(jfp)
c              end do
           endif
            jj=1
            do jfp=1,nfp
               call vdot(pd(jj),1,pd(jj),1,tpd(jfp),nmp)
               jj=jj+nmpfp
            end do
            if (debug) write(*,122)'tpd=',(tpd(jtg),jtg=n1,n4),
     &      (tpd(jtg),jtg=1,ntg)
c  Solve matrix equation PD.dx=resids
           call vmov(zero,0,wk,1,nfp)
           call shfti(pd,nmpfp,nmpfp,nfp,resids,nmpfp,
     &     i1,tau,krank,rdum,wk,ip)
c
           call vmov(resids,1,dx,1,nfp)
           if(debug) then
           write(*,122)'dx=',(dx(jtg),jtg=n1,n4),(dx(jtg),jtg=1,ntg)
              if(krank.lt.nfp) then
                  write(6,*)'Rank Deficient:',krank,'  /',nfp
              else
                  write(6,*)'Full Rank:',krank,'  /',nfp
              endif
           endif
           if(dx(n1).gt.9.21) then
             fs=10000.
           elseif(dx(n1).lt.-9.21) then
             fs=.0001
           else
             fs=exp(dx(n1))
           endif
           dx(n1)=cx(n1)*(fs-1)
        endif  ! (rms.gt.9*rwas)
        if(kconv.eq.0 .or. mit.eq.0) go to 63   ! convergence
c  Limit maximum step size for fitted gases (non-linear).
c        do jtg=1,ntg
c        dxlimit=0.2+abs(cx(jtg))
c        if(abs(dx(jtg)) .gt. dxlimit) dx(jtg)=sign(dxlimit,dx(jtg))
c        end do
c        write(*,*)'nit,dx=',nit,dx
c        call vadd(cx,1,dx,1,cx,1,nfp)
        fr=1.0
        do jfp=1,nfp
           dxlimit=0.2+abs(cx(jfp))
           xfr=abs(dxlimit/abs(dx(jfp)))
           if(xfr.lt.fr) then
              fr=xfr
              kfp=jfp
           endif
        end do
        if(debug .and. fr.lt.1.)
     &  write(*,*)'Limited step size: kfp,fr=',kfp,fr
        call vsma(dx,1,fr,cx,1,cx,1,nfp)
        rwas=rms
        shift=gint*dble(cx(n3)-ax(n3))
 65     continue
 63     totit=totit+nit
        cx(n3)=ax(n3)+sngl(shift/gint)
c=========================================================================
c  Compute upper-triangle of covariance matrix & move diagonal elements into EX
       var=(rms*corrld)**2/(1.-float(nfp)/(float(nmp)+0.1))
       var=var+tau**2  ! fudge to prevent EX=0 when rms=0
       call scov2(pd,nmpfp,nfp,ip,var,ierr)
       if(debug) then
       write(*,*)' Upper triangle of correlation matrix:'
       write(*,'(8(a,i1))')('  Target_',j,j=1,ntg),
     & '   Cntuum     Tilt      Shift     ZOff'
       do i=1,nfp
          write(*,'(11f10.6)')(pd((j-1)*nmpfp+i)
     & /sqrt(pd((j-1)*nmpfp+j)*pd((i-1)*nmpfp+i)),j=i,nfp)
       end do
       endif
       if(ierr.eq.0 .and. krank.gt.0) then
          call vmov(pd,nmpfp+1,ex,i1,nfp)
       else
          call vmov(1/tau**2,0,ex,i1,nfp)
       endif
       ex(n1)=ex(n1)*cx(n1)**2
c       write(6,*)(cx(j),j=1,nfp)
c       write(6,*)(dx(j),j=1,nfp)
c       write(6,*)(sqrt(ex(j)),j=1,nfp)
c=========================================================================
c  Add RSS error contributions from change in zero level (ZERR) and RMS fit
c       call vmul(dx,1,dx,1,dx,1,nfp)
c       call vadd (dx,1,ex,1,ex,1,nfp)
c       call vadd (ex,1,(3*rms/cx(n1))**2,0,ex,1,nfp)   !  This is a fudge
c       call vsqrt(ex,1,ex,1,nfp)
       cl=cx(n1)
       if(abs(cl).lt.tiny) cl=tiny
       do jfp=1,nfp
          ex(jfp)=sqrt(abs(ex(jfp))
c     &    +dx(jfp)**2      ! perturbation from adding ERROFF to residuals
     &    +(3*rms/cl)**2   ! fudge
c     &    +(100*(cx(jfp)-ax(jfp))*rms/cl)**2   ! fudge
c     &    +25*dxwas(jfp)**2   ! fudge
     &    +25*dx(jfp)**2   ! fudge
c     &    +0.04*(cl**2+1./cl**2)*(cx(jfp)-ax(jfp))**4  ! fudge
     &    )
       end do
c===================================================================
c  Compute transmittances (convolved with ILS) of individual target gases.
c  Place them in PD, which just so happens to be exactly the right size.
96    jva=1
      jpd=1
      kn3=1+ncp*(n2) ! Start address of workspace
      sh=rdec*cx(n3)
      do 95 jtg=0,ntg
        if(jtg.eq.0) then ! non-target gases
          call vexp(spxv(jva),i1,spxv(kn3),i1,ncp)
        else        ! target gases
          call vmul(spxv(jva),i1,cx(jtg),i0,spxv(kn3),i1,ncp)
          call vexp(spxv(kn3),i1,spxv(kn3),i1,ncp)
        endif
        call newdec(spxv(kn3),ncp,slit,nii,ldec,rdec,sh,pd(jpd),nmp)
        jva=jva+ncp
        jpd=jpd+nmpfp
95    continue
c
c  Do the solar spectrum too. Put them in pd(*,n2) for convenience.
      call newdec(spxv(jva),ncp,slit,nii,ldec,rdec,sh,pd(jpd),nmp)
c===========================================================================
c  The following obscure section of code computes the vertical column abundance.
c  It should really be removed from GFIT since
c  the vertical column cannot be uniquely determined from the fitted spectrum.
      if (nlev.le.2) then
        spver(1)=splos(1)
        spver(nlev)=splos(nlev)
      else
        call vmov(zero,0,spver,1,mlev)
c        do klev=2,nlev  !   Changed Mar 12, 2006  GCT
        do klev=3,nlev
           if(z(klev).gt.zmin) go to 777
        end do
        write(6,*) 'Warning: zmin exceeds z(nlev)',zmin,z(nlev)
 777    dz=z(klev)-z(klev-1)
        xo=(zmin-z(klev-1))/dz
        if(d(klev).le.0.0) then
           logrp=0.0
        else
           logrp=log(d(klev-1)/d(klev))
        endif
c        write(6,*)'xo=',klev,xo,zmin,z(klev),z(klev-1),d(klev-1),d(klev)
        xl=logrp*(1.-xo)
        spver(klev-1)=dz*(1.-xo)*
     &  (1-xo-xl*(1+2*xo)/3+xl**2*(1+3*xo)/12+xl**3*(1+4*xo)/60)/2
        spver(klev)  =dz*(1.-xo)*
     &  (1+xo+xl*(1+2*xo)/3+xl**2*(1+3*xo)/12-xl**3*(1+4*xo)/60)/2
        do jlev=klev+1,nlev
          dz=z(jlev)-z(jlev-1)
          if(d(jlev).le.0.0) stop 'Non-positive density'
          logrp=log(d(jlev-1)/d(jlev))
          spver(jlev-1)=spver(jlev-1)+
     &    dz*(1.-logrp/3+logrp**2/12-logrp**3/60)/2
          spver(jlev)  = dz*(1.+logrp/3+logrp**2/12+logrp**3/60)/2
        end do
      endif  ! nlev.le.2
c===================================================================
c Compute vertical and LOS slant columns for the target gases
c Do in reverse order so that the concentration profile of the
c first target gas is retained in CP at the end.
      do jtg=ntg,1,-1
         jspexi=speci(jtg)
         call vmul(vmr(jspexi),nspexi,d,i1,cp,i1,nlev)
         call vdot(cp,i1,spver,i1,overcol(jtg),nlev)
         call vdot(cp,i1,splos,i1,oloscol(jtg),nlev)
         airmass(jtg)=oloscol(jtg)/overcol(jtg)
         if(airmass(jtg).gt.999.999) airmass(jtg)=999.999
         if(airmass(jtg).lt.-99.999) airmass(jtg)=-99.999
         overcol(jtg)=overcol(jtg)*ckm2cm
         oloscol(jtg)=oloscol(jtg)*ckm2cm
      end do
c
      if(ntg.eq.0) then
c   compute column amounts for air (vmr=1.000)
         call vdot(d,1,spver,1,overcol(1),nlev)
         call vdot(d,1,splos,1,oloscol(1),nlev)
         airmass(1)=oloscol(1)/overcol(1)
         overcol(jtg)=overcol(jtg)*ckm2cm
         oloscol(jtg)=oloscol(jtg)*ckm2cm
c         write(29,76)' '//runlab,0,0.0,0.0,0.0,
c     &   0.0,0.0,zmin,tairmass, tvercol,1.0,0.0
       endif
c
       call vdot(curv,1,resids,1,cr,nmp)
       call vdot(curv,1,curv,1,cc,nmp)
       cr=0.05*cr*(1+amal/fovi)**3/cc/frqcen*(gint/rect)**2
       if(cr.ge.amal**2) cr=amal**2
       misalign=1000*sqrt(abs(amal**2-cr))/100
c
c This section of code added by GCT 050113 to calculate solar shift
c Compute solar and gas shifts. Modified 20060106.
       tns=0.0
       tds=tiny
       tng=0.0
       tdg=tiny
       if(index(winfo,' so ').gt.0) then
       nop=(nmp+1)/2
       do k=2,nmp-1
          cntuum=cx(n1)*(1.+cx(n2)*float(k-nop)/(nmp-1))
          residual=obsrvd(k)-calcul(k)
c
c Solar Shift (pd(k+n1*nmpfp) is the solar transmittance spectrum)
          gs=cntuum*(pd(k+1+n1*nmpfp)-pd(k-1+n1*nmpfp))
          tns=tns+gs*residual
          tds=tds+gs**2
c
c Gas shift  (calcul/solar is the gaseous transmittance spectrum)
          gg=cntuum*
     &(calcul(k+1)/pd(k+1+n1*nmpfp)-calcul(k-1)/pd(k-1+n1*nmpfp))
          tng=tng+gg*residual
          tdg=tdg+gg**2
       end do     ! k=2,nmp-1
c       write(*,*) 'solar-gas shift=',tns,tds,tng,tdg,residual,cntuum,
c     & 2*gint*(tns/tds-tng/tdg)/frqcen
       endif

c  Solar-Gas shift is the difference of the solar shift and the gas shift.
c  Divide by frequency to normalize to a dimensionless stretch which
c  should be the same for all windows. Multiplied by 10^6 for ppm.
       sgshift=1E6*2*gint*(tns/tds-tng/tdg)/frqcen

       wsgshift=sgshift
       if(sgshift.gt.+99.9) wsgshift=99.9
       if(sgshift.lt.-9.9) wsgshift=-9.9

       avgrms=avgrms+abs(rms)
       avgcl=avgcl+abs(cx(n1))
       toterr=toterr+1.0/ex(1)**2
       tottc=tottc+cx(1)/ex(1)**2
       tottc2=tottc2+((cx(1)-1)/ex(1))**2

       wrms=abs(100.*rms/cx(n1))
       if(wrms.gt.+9.9999) wrms=9.9999

       wshift=1000*shift
       if(wshift.gt.+99.9) wshift=99.9
       if(wshift.lt.-9.9) wshift=-9.9

       wtilt=100*cx(n2)
       if(wtilt.gt.+99.9) wtilt=99.9
       if(wtilt.lt.-9.9) wtilt=-9.9

       wcl=cx(n1)
       if(wcl.gt.+9.999) wcl=9.999
       if(wcl.lt.-.999) wcl=-.999

       wcle=abs(ex(n1))
       if(wcle.gt..999) wcle=.999

       wzo=cx(n4)
       if(wzo.gt.+9.999) wzo=9.999
       if(wzo.lt.-.999) wzo=-.999

       if(zmin.gt.999.999)zmin=999.999

       do jtg=1,ntg
         if(cx(jtg).gt.9999.9999) cx(jtg)=9999.999
         if(cx(jtg).lt.-999.9999) cx(jtg)=-999.9999
       end do
       write(6,76)runlab,nit,wcl,
     & wtilt,wshift,wsgshift,wzo,wrms,zmin,(airmass(jtg),
     & overcol(jtg),cx(jtg),ex(jtg),jtg=1,min0(1,ntg))
       write(lunt,76)runlab,nit,wcl,
     & wtilt,wshift,wsgshift,wzo,wrms,zmin,(airmass(jtg),
     & overcol(jtg),cx(jtg),ex(jtg),jtg=1,ntg)
 76    format(1x,a21,i2,1x,f5.3,3(1x,f4.1),1x,f5.3,
     & 1x,f6.4,f7.3,9(0pf7.3,1pe10.3,0pf9.4,1pe8.1))
c
c  Write .spt file for the first MSPT spectral fits
       if(ispec .lt. mspt) then
          call write_spt(19,winfo,sptfile(:lnbc(sptfile))//runlab,
     &    obsrvd,calcul,cx,ex, startm+shift, dopp*freq_flag,
     &    gint,overcol, pars, asza+zenoff,obalt,zmin,wrms,
     &    pd, nmp, nmpfp, ntg)
c
c  Output PD's (weighting functions), for subsequent use
c  in deriving averaging kernels.
        if( index(winfo,' ak ') .gt. 0) then
           open(luna,file=akfile(:lnbc(akfile))//runlab,
     &     status='unknown')
           call ak(luna,slit,nii,ldec,spxv,splos,vac,nlev,ncp,rdec,
     &      cx,ntg,pd,nmp)
c  Level 1 is the cell, so start jlev at 2.
           write(luna,*)(splos(jlev)*cp(jlev),jlev=2,nlev)
           write(luna,*) pout/1013.25
           write(luna,*)(p(jlev),jlev=2,nlev)
           close(luna)
        endif
       endif ! (ispec.lt.2500)
      end do   !  ispec=1,nspec     Main fitting loop over spectra
      write(6,*)' Total number of iterations =',totit
      write(6,*)' Total number of spectra fitted  =',nspec
      write(6,*)' Total number of mavblocks found =',mavfound
      write(6,*)' Number of spectra w/o mavblocks =',mavwithout
      write(6,*)' Average % RMS fit =',100*avgrms/avgcl
      avgtc=tottc/toterr
      rmstc=sqrt(abs(tottc2/toterr-(avgtc-1.)**2))
      write(6,*)' Average Totcon =',avgtc,' +- ',rmstc
      close(lunq)
      close(lunr)
      close(lunt)
      stop
      end

      function riair(w,t,p,h)
c  Edlen's formula for the refractive index of air at wavenumber w cm-1
c  T is the air temperature in degrees Celsius
c  P is the air pressure in mbar
C  H is the relative humidity in %
      real*8 w,t,p,h,pp,delt,hh,f1f,f2f
      real*8 riair
      F1F(w)=0.378125+w*w*(2.1414E-11+w*w*1.793E-21)
      F2F(w)=0.0624-w*w*6.8E-14
      if(w.lt.0.0) w=0.0
      PP=(p-.3175+5.E-4*(p-745.)+.13*(t-10.))*.7500646
      DELT=PP*(1.-PP*(1.57E-8*t-1.049E-6))/(1.+3.661E-3*t)
      HH=h*EXP(1.52334+t*(.07217-2.9549E-4*t))/(100.+.3661*t)
      riair=1.D0+1.E-6*(DELT*F1F(w)-HH*F2F(w))
      return
      END
