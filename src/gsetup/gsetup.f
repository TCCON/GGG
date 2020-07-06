c   gsetup.f
c   Provides user interface to create input files for gfit
c
c  Inputs:
c     Menu-driven user input
c
c  Outputs:
c     multigg.sh           ! Batch file
c     gas_1234_runlog.ggg  ! GFIT input files (one for each window)
c     runlog.mav           ! P/T/vmr profiles interpolated onto chosen levels
c     runlog.ray           ! Ray-traced atmospheric slant paths
c     gsetup.rpt           ! File containing report
c
c  
c
c  The following pseudo-code summarizes the main program structure
c      nhr=3 or 24 (model time-step in hours)
c      modname_bak='jplsummr.mod                                    '
c      vmrname_bak='summer_35N.vmr                                  '
c      Choose geometry,runlog,levels,windows
c      do kwin=1,mwin
c         write .ggg file for each window
c      end do
c      modname_cur='xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
c      read_runlog_header
c      do ispec = 0,999999
c         read_runlog:  iyr,idoy,zpdtim,oblat,oblon
c         jul = julian(iyr,1,idoy)
c         iyyy,imm,idd = caldat(jul)
c         ihh = nhr*nint(zpdtim/nhr)
c         if(ihh.lt.0) then
c            ihh=ihh+24
c            idd=idd-1
c         endif
c         if(ihh.ge.24) then
c            ihh=ihh-24
c            idd=idd+1
c         endif
c
c         modname_new = ‘FPIT_’,iyyy,imm,idd,ihh,
c      &  nint(abs(lat)),ns,nint(abs(lon)),ew,’.mod’
c         newmod=.false.  ! Current model okay
c         if(modname_new.ne.modname_cur) then
c            newmod_exist=inquire(modname_new)
c            if(newmod_exist .eqv. .true.) then
c               modname_cur=modname_new
c               nmod=nmod+1
c            else
c               if(fpit) stop 'Missing .mod file'
c               modname_cur=modname_bak
c            endif
c            newmod=.true.
c            call read_mod_fc(modname_cur)
c         endif  !       if(modname_new.ne.modname_cur) then
c 
c         vmrname_new=specname(1:10)//'.vmr' !  TCCON/NCEP
c         vmrname_new='JL1_YYYYMMDDHH'//'.vmr' ! FPIT
c         newvmr=.false.
c         if(vmrname_new.ne.vmrname_cur) then
c            newvmr_exist=inquire(vmrname_new)
c            if(newvmr_exist .eqv. .true.) then
c               vmrname_cur=vmrname_new
c               nvmr=nvmr+1
c            else
c               if(fpit) stop 'Missing vmr file'
c               vmrname_cur=vmrname_bak    ! NCEP
c            endif
c            newvmr=.true.
c            call read_refvmrs(vmrname_cur)
c            if ((nlev.gt.1).and.(modtype.ne.'FPIT')) then
c               apvmr=adjust(refvmrs) ! old code for NCEP
c            else
c               apvmr=refvmr  ! no adjustment for FPIT case
c            endif
c         endif     !   if(vmrname_new.ne.vmrname_cur)
c
c         if(newmod.or.newvmr) then
c            write(.mav file) specname(ispec), z,t,p,d,vmr()
c            call write_mav(apvmr)
c         endif
c         call tlpath(…….splos)
c         write(ray_file) specname,…splos(j),j=1,nlev)
c
c      end do  !      do  ispec = 0,999999

c
c  There are 3 modname_xxx and vmrname_xxx suffixes
c      name_bak ! the back-up mod/vmr file name (summer_35N)
c      name_new ! the desired new mod/vmr file.
c      name_cur ! the current mod/vmr file name (_new or _bak)

      implicit none
      real*4 spi,d2r          !PI in single precision
      parameter(spi=3.14159265,d2r=spi/180.)

      character md5sum*24
      parameter(md5sum = "$GGGPATH/bin/gfit_md5sum")

      include "../gfit/ggg_int_params.f"

      integer*4 bytepw,fbc,fnbc,lnbc,ifirst,ilast,ilev,
     & idoy,iyr,j,w1,w2,le,lg,lc,lrt,idum,
     & nlhead,ncol_iso,ncol_win,mwin,kwin,jwin,kk,
c     & ih,ihp,ihz,
     & lunw_ray,lunw_ggg,
     & lunr_rlg,lunw_mul,
     & lunw_mav,lunw_rpt,lunr_mod,
     & lunr_vmr,
     & lunr_csi,
     & lunr_iso,
     & lunr_lev,
     & lunr_win,
     & jgas,sl,
     & jd,
     & nlev,nlhead_ray,nlhead_ggg,nlhead_win,nspeci,ngas,
     & lr,lspn,i,doy,iyyyy,imm,idd,ihh,jhh,jul,jwas,nhr,
     & multipath,nmod,nvmr,ninqmod,ninqvmr,jspe,possp,rc,istat,lw,
c     & lm,lv,
     & lunr_men, lunr_stdin,system
      parameter (lunw_ggg=56)    ! for writing the .ggg files
      parameter (lunr_iso=57)    ! for reading isotopologs.dat
      parameter (lunw_mav=60)    ! for writing the .mav file
      parameter (lunr_men=61)    ! for reading the various .men files
      parameter (lunr_lev=62)    ! for reading the levels files
      parameter (lunr_mod=63)    ! for reading the .mod file (readmodFC)
      parameter (lunr_win=64)    ! for reading the windows file
      parameter (lunw_mul=66)    ! for writing the multiggg.sh file
      parameter (lunw_ray=68)    ! for writing the .ray file
      parameter (lunr_rlg=69)    ! for reading the runlog
      parameter (lunw_rpt=72)    ! for writing gsetup.rpt
      parameter (lunr_vmr=73)    ! for reading the .vmr file (readvmrFC)
      parameter (lunr_csi=74)    ! for reading xx_cell_status_info.dat
      parameter (nlhead_ray=4)   ! Number of header lines in .ray file.
      parameter (nlhead_ggg=18)  ! Number of header lines in .ggg file.
      parameter (mwin=800)       ! Max number of windows

      logical*4 newmod,newvmr,newmod_exist,newvmr_exist

      integer*4
     & igas_in_cell(ncell),kcell

      real*4 dum,frqcen,height,tlat,compute_ztrop,
     & compute_ztrop_gct,ztrop_gct,delta_p,
     & zmin,zpres,
c     & histop(100),dhp,
c     & histoz(100),dhz,
     & p_cell(ncell),t_cell(ncell),vmr_cell(ncell),
     & zobs,               ! the geocentric radius of the observer fed to TLPATH
     & bend,               ! total ray bending due to refraction
     & d(mlev),            ! number densities at levels (molec.cm-3)
     & h2o_dmf(mlev),      ! h2o profile from model file
c     & co2vmr(mlev),      ! co2 profile from simulate_co2_vmr
     & fovr,               ! External angular radius of FOV in degrees
     & mmw(mlev),          ! mean molecular weights
     & p(mlev),            ! pressures of levels (atm.)
c     & r(mlev),            ! radii of levels (km) = z + roc
     & roc,                ! radius of curvature (km)
     & rocx,               ! radius of curvature (km)
     & cell_length(ncell), ! Length of cell permanently in solar beam
     & solzen,             ! solar zenith angle used (may be refracted)
     & splos(mlev),        ! array of line-of-sight slant paths
     & t(mlev),            ! temperatures of levels (K)
     & gradlat(mgas),      ! Latitude gradients of each gas
     & seacycle(mgas),     ! Seasonal cycle of each gas
     & strend(mgas),       ! Secular Trend of each gas
     & refvmr(mgas,mlev),  ! buffer for reference vmr's
     & apvmr(mgas,mlev),   ! buffer for a priori vmr's
     & z(mlev)             ! altitudes of levels (km)

      real*8  fryr,latwas,lonwas,diff,
     & itcz_lat,itcz_width,
     & reflat_vmr,date_mod,date_vmr,ztrop_vmr,ztrop_mod

      real*8  zpbl, wlimit, gcd,
     & trlg,
     & frdoy,frwas,        ! Fractional Julian day numbers
c     & eta, eta2,
     & oblat,              ! observation latitude (deg).
c     & coblat,             ! observation latitude (deg).
c     & soblat,             ! observation latitude (deg).
c     & coblat2,            ! observation latitude (deg).
c     & soblat2,            ! observation latitude (deg).
     & oblon,              ! observation longitude (deg).
     & obalt,              ! observation altitude (km)
c     & tplat,              ! observation latitude (deg).
c     & tplon,              ! observation longitude (deg).
c     & tpalt,              ! observation altitude (km)
     & asza,               ! astronomical solar zenith angle (unrefracted)
c     & asza1,              ! astronomical solar zenith angle (unrefracted)
c     & eorv,ervc,
     & opd,                ! Optical path difference (cm) of interferogram
     & graw,               ! spacing of raw spectrum (cm-1) from GETINFO
     & zpdtim, zpdwas,     ! Time of ZPD (UT hours)
     & zenoff,             ! Zenith angle pointing offset (deg)
     & azim,               ! Solar Azimuth Angle
c     & azim1,              ! Solar Azimuth Angle
     & osds,               ! Observer-Sun Doppler Stretch (ppm_
     & fovi,               ! Internal angular diameter of FOV (radians)
     & fovo,               ! External angular diameter of FOV (radians)
     & amal,               ! angular misalignment of interferometer (radians)
     & zoff,               ! Zero level offset (dimensionless fraction)
     & snr,                ! Signal-to-Noise Ratio (dimensionless)
     & tins,               ! Inside temperature
     & pins,               ! Inside pressure
     & hins,               ! Inside humidity
     & tout,               ! Outside temperature
     & pout,               ! Outside pressure
     & hout,               ! Outside humidity
     & sia,                ! Solar Intensity Average (arbitrary units)
     & fvsi,               ! Fractional Variation in Solar Intensity
     & wspd,               ! Wind Speed (m/s)
     & wdir,               ! Wind Direction (deg)
     & aipl,               ! Airmass-Independent Path length (km)
     & lasf,               ! Laser Frequency (e.g. 15798 cm-1)
     & wavtkr              ! Frequency of the sun-tracker (actively-tracked)

        
      character
     & window_str*160,       ! window string
     & filnam*128,           ! general purpose character buffer
     & header*92,            ! general purpose character buffer
     & data_fmt_read_rl*256,
     & col_labels_rl*320,
     & levels*64,            ! name of file containing levels
     & listof*64,            ! the selected list of windows
     & newmkivname*10,       ! del
     & modname_bak*80,       ! name of backup model
     & modname_cur*80,       ! name of current model
     & modname_new*80,       ! name of desired new model
     & vmrname_bak*80,       ! name of backup vmr file
     & vmrname_cur*80,       ! name of current vmr file
     & vmrname_new*80,       ! name of new vmr file
     & menupath*(mfilepath), ! path to xxx.men file
     & prvwin*14,            ! name of previous window (i.e. GAS_1234)
     & llfile*128,           ! Path to a linelist
c    & user*8,               ! investigator
     & ray_header_fmt*26,    ! Format statement in .ray file header
     & isofile*(mfilepath),  ! path to "isotopolog.dat"
     & vmrlabel*1024,        ! column labels from vmr file
     & str_isotop*294        ! line from isotopologs.dat file

      character
     & modtype*4,
     & tll_file*80,          ! Telluric Linelists
     & winnam(mwin)*32,      ! Windows (gas_1234)
     & col1*1,               ! first column of runlog record
     & apf*2,                ! apodization function (e.g. BX N2, etc)
     & dl*1,                 ! forward or backward slash
     & ext*3,                ! geometry ['air','bal','gnd','lab',orb','syn']
     & ns*1, ew*1,
     & gggdir*(mpath),       ! ggg directory path (GGGPATH?)
     & specname*(nchar),     ! spectrum name
     & version*64,           ! current program version
     & ray_data_fmt*38,      ! format of data in .ray file
     & rlgfile*120,          ! name of runlog file
     & menupathutfile*20     ! gsetup.input

c      data histop/100*0.0/
c      data histoz/100*0.0/
c      dhp=10.  ! pressure bin width (mbar)
c      dhz=0.2  ! altitude bin width (km)

      data strend/mgas*0.0/
      data gradlat/mgas*0.0/

      version=
     & ' GSETUP                   Version 4.70        2020-06-29   GCT '

      ninqmod=0
      ninqvmr=0
      lspn=0      ! avoid compiler warning (may be used uninitialized)
      modname_bak='jplsummr.mod                                    '
      vmrname_bak='summer_35N.vmr                                  '
      ray_data_fmt='(a57,6(1x,f11.5),2f12.6,260(1x,f11.5))'
      zpdtim=-9999.9
      oblat=-9999.9
      oblon=-9999.9
      frdoy=-9999.9
      frwas=frdoy  ! avoid compiler warning (may be used uninitialized)
      latwas=oblat ! avoid compiler warning (may be used uninitialized)
      lonwas=oblon ! avoid compiler warning (may be used uninitialized)
      zpdwas=zpdtim

      idum=mauxcol ! prevent compiler warning (unused variable)
      idum=mcolvav ! prevent compiler warning (unused variable)
      idum=mrow_qc ! prevent compiler warning (unused variable)
      idum=mspeci  ! prevent compiler warning (unused variable)
      idum=mvmode  ! prevent compiler warning (unused variable)

      jwas=-1

c     Platform specification:      DG090519
      call get_ggg_environment(gggdir, dl)
      lrt=lnbc(gggdir)     !Length of gggdir
c
c  Choose which telluric linelists to be used
      tll_file='telluric_linelists.md5'
      llfile=gggdir(:lrt)//'linelist'//dl//'atm.161'
      istat=system(md5sum//' '//llfile//' > '//tll_file)
      llfile=gggdir(:lrt)//'linelist'//dl//'atmnv.275'
      istat=system(md5sum//' '//llfile//' >> '//tll_file)
      llfile= gggdir(:lrt)//'linelist'//dl//'pll.101'
      istat=system(md5sum//' '//llfile//' >> '//tll_file)
      llfile=gggdir(:lrt)//'linelist'//dl//'fcia.101'
      istat=system(md5sum//' '//llfile//' >> '//tll_file)
      llfile=gggdir(:lrt)//'linelist'//dl//'scia.101'
      istat=system(md5sum//' '//llfile//' >> '//tll_file)

c     llfile=gggdir(:lrt)//'linelist'//dl//'hitran00.101'
c     llfile=gggdir(:lrt)//'linelist'//dl//'hitran04.161'
c     llfile=gggdir(:lrt)//'linelist'//dl//'hitran08.161'
c     llfile=gggdir(:lrt)//'linelist'//dl//
c     & 'ab-initio-h2o-8000_hitran.161'
c     llfile=gggdir(:lrt)//'linelist'//dl//'H2O-GEISA-hitraned.101'
c     llfile=gggdir(:lrt)//'linelist'//dl//'16MiKaMo_hitran.161 '
c     llfile=gggdir(:lrt)//'linelist'//dl//'stripped03_hitran08.161 '
c     llfile=gggdir(:lrt)//'linelist'//dl//'stripped03_hitran13.161 '
c     llfile=gggdir(:lrt)//'linelist'//dl//'03_hitran16_beta.162 '
c     llfile=gggdir(:lrt)//'linelist'//dl//'stripped03_atm.101 '
c     llfile=gggdir(:lrt)//'linelist'//dl//'leftover03_atm.101 '
c     llfile=gggdir(:lrt)//'linelist'//dl//'stripped12_hitran13.161 '
c     llfile=gggdir(:lrt)//'linelist'//dl//'leftover12_atm.161 '
c
c  Find NSPECI
      isofile=gggdir(:lrt)//'isotopologs/isotopologs.dat'
c      isofile='isotopologs_local.dat'
      open(lunr_iso,file=isofile,status='old')
      read(lunr_iso,*) nlhead,ncol_iso
      do i=2,999
1        read(lunr_iso,'(a)',end=76) str_isotop
         if(str_isotop(1:1).eq.';' .or. str_isotop(1:1).eq.':') go to 1
      end do
76    nspeci=i-1-nlhead
      close(lunr_iso)

      ext(1:3)='   '
      open(lunw_rpt,file='gsetup.rpt', status='unknown')

c      write(62,*) 8,4
c      write(62,*) 'fort.62  compares Ztrop computed on Pgrid and Zgrid'
c      write(62,*) 'fort.62  compares Ztrop computed on Pgrid and Zgrid'
c      write(62,*) 'fort.62  compares Ztrop computed on Pgrid and Zgrid'
c      write(62,*) 'fort.62  compares Ztrop computed on Pgrid and Zgrid'
c      write(62,*) 'fort.62  compares Ztrop computed on Pgrid and Zgrid'
c      write(62,*) 'fort.62  compares Ztrop computed on Pgrid and Zgrid'
c      write(62,*) ' date  year  ztrop_pgrid  ztrop_zgrid'

c      write(83,*)3,6
c      write(83,*)' fort.83 summarizes tropopause information'
c      write(83,*)'date ptrop_ncep ztrop_ncep ztrop_ncep2 ptrop_gct
c     &  ztrop_mod'

c      write(84,*)3,3
c      write(84,*)' fort.84 is a histogram of NCEP tropopause pressure'
c      write(84,*)' i  pressure  number_of_days'

c      write(85,*)3,3
c      write(85,*)' fort.85 is a histogram of NCEP tropopause altitude'
c      write(85,*)' i  pressure  number_of_days'

      write(6,*) version

c  choose an observation geometry
      if (iargc() == 0) then
         write(6,9913)
 9913    format(' Geometry (a=air,b=bal,g=gnd,l=lab,o=orb,s=syn) ? ',$)
         read(5,'(a)') ext(1:1)
         lunr_stdin = 5
      elseif (iargc() == 1) then
         call getarg(1, menupathutfile)
         lunr_stdin = 10
         open(lunr_stdin, file=menupathutfile, status='old')
         read(lunr_stdin,'(a)') ext(1:1)
      else
         write(*,*) 'Usage: $gggpath/bin/gsetup inputfile_containing_'//
     &   'selections'
         stop
      endif
      if(ext(1:1).eq.'a') ext(1:3)='air'
      if(ext(1:1).eq.'b') ext(1:3)='bal'
      if(ext(1:1).eq.'g') ext(1:3)='gnd'
      if(ext(1:1).eq.'l') ext(1:3)='lab'
      if(ext(1:1).eq.'o') ext(1:3)='orb'
      if(ext(1:1).eq.'s') ext(1:3)='syn'
      if(ext(2:3).eq.'  ') stop 'unknown geometry'

c------------------------------------------------------------------
c  Choose which runlog to analyze
      menupath=gggdir(:lrt)//'runlogs'//dl//ext//dl//'runlogs.men'
      call readmenu(lunr_men,menupath,lunr_stdin,rlgfile)
      lr=lnbc(rlgfile)

c------------------------------------------------------------------
c  choose which levels to use
      menupath=gggdir(:lrt)//'levels'//dl//'levels.men'
      call readmenu(lunr_men,menupath,lunr_stdin,levels)
c  Read chosen level altitudes
      open(lunr_lev,file=gggdir(:lrt)//'levels'//dl//
     & levels(:lnbc(levels)),status='old')
      do ilev=1,mlev
         read(lunr_lev,*,end=877) z(ilev),mmw(ilev)
      end do
      read(lunr_lev,*,end=877) dum
      write(6,*) 'Warning: number of levels exceeds NLMAX'
 877  nlev=ilev-1
      close(lunr_lev)
      if(nlev+1.gt.mlev) stop 'Error: nlev+1 > mlev '
c------------------------------------------------------------------
c  choose which list of windows to analyze
      menupath=gggdir(:lrt)//'windows'//dl//ext//dl//'windows.men'
      call readmenu(lunr_men,menupath,lunr_stdin,listof)

c  Choose model type
      write(6,9914)
 9914 format(' Std TCCON processing: 3-hourly model & VMR (y/n) ?',$)
      read(lunr_stdin,'(a)')modtype

      if(modtype.eq.'y' .or. modtype.eq.'Y' ) then
         if (levels(:lnbc(levels)).ne.'ap_51_level_0_to_70km.gnd') stop
     &   ' ap_51_...  levels file must be used for TCCON processing'
         nhr=3
         modtype='FPIT'
      else
         nhr=24
         modtype='NCEP'
      endif

      close(lunr_stdin)

c  read the individual microwindows and create the xxxxxxxx.ggg input files
      if(dl.eq.char(92))then  ! Back-Slash (\)
         open(lunw_mul,file='multiggg.bat',status='unknown')
      else
         open(lunw_mul,file='multiggg.sh',status='unknown')
      endif
      open(lunr_win,file=gggdir(:lrt)//'windows'//dl//ext
     $ //dl//listof(:lnbc(listof)),status='old')
      read(lunr_win,*)nlhead_win,ncol_win
      do j=2,nlhead_win
         read(lunr_win,'(a91)') header
      end do
      prvwin='            '
      do kwin=1,mwin
 666     read(lunr_win,'(a)',end=103) window_str
         if(window_str(1:1).eq.':' .or. window_str(1:1).eq.';') goto 666
         w1=fnbc(window_str)
         if(w1.eq.0) go to 666 ! blank line
         w2=index(window_str,'.')
c         write(*,*) w1,w2,window_str(w1:w2)

         lc=index(window_str,':')  ! lc=position of colon
         lg=lc+fnbc(window_str(lc+1:)) ! start of first target gas
         le=lg+fbc(window_str(lg+1:))-1  ! end of first target gas
         if(le.le.lc) then  ! No target gas
            winnam(kwin)='__'//window_str(w1:w2)
         else
            winnam(kwin)=window_str(lg:le)//'_'//window_str(w1:w2)
         endif
         lw=lnbc(winnam(kwin))
         filnam=winnam(kwin)(:lw)//rlgfile(:lr-3)//'ggg'

         kk=0   ! Number of previous, identically-named windows
         do jwin=1,kwin-1
            if(winnam(jwin).eq.winnam(kwin)) kk=kk+1
         end do
         if(kk.gt.0) filnam=filnam(:lw-1)//char(96+kk)//filnam(lw:)

c         write(*,'(3i3,a)') lc,lg,le,' '//window_str(:lnbc(window_str))
         open(lunw_ggg,file=filnam,status='unknown')
         write(lunw_ggg,'(i3)') nlhead_ggg
         write(lunw_ggg,'(a)') version
         write(lunw_ggg,'(a)')
     &       gggdir(:lrt)//'config'//dl//'data_part.lst'
         write(lunw_ggg,'(a)')
     &       gggdir(:lrt)//'apriori'//dl//'gfit_ap.'//ext
         write(lunw_ggg,'(a)') gggdir(:lrt)//'runlogs'//dl//ext//dl//
     $   rlgfile(:lr)
         write(lunw_ggg,'(a)') gggdir(:lrt)//'levels'//dl//levels
         write(lunw_ggg,'(a)') gggdir(:lrt)//'models'//dl//ext//dl
         write(lunw_ggg,'(a)') gggdir(:lrt)//'vmrs'//dl//ext//dl
         write(lunw_ggg,'(a)') rlgfile(:lr-3)//'mav'
         write(lunw_ggg,'(a)') rlgfile(:lr-3)//'ray'
         write(lunw_ggg,'(a)') isofile(:lnbc(isofile))
         write(lunw_ggg,'(a)')
     &   gggdir(:lrt)//'windows'//dl//ext//dl//listof

c  write names of linelists

c         write(lunw_ggg,'(a)')      ! llsize
         write(lunw_ggg,'(a)') tll_file
         if(index(window_str(w2:),' sg ').gt.0 .or.
     &   index(window_str(w2:),' so/').gt.0   .or.
     &   index(window_str(w2:),' so ').gt.0 ) then
            write(lunw_ggg,'(a)')
     &      gggdir(:lrt)//'linelist'//dl//'solar_merged.108'
         else
            write(lunw_ggg,'(a)')
         endif

c  write location of $gggpath/ak/jxxxxx files (Jacobians)
         write(lunw_ggg,'(a)') gggdir(:lrt)//'ak'//dl//'j'

c  write location of spectral fits
         write(lunw_ggg,'(a)') gggdir(:lrt)//'spt'//dl//'z'

         write(lunw_ggg,'(a)') filnam(:lnbc(filnam)-3)//'col'
         write(lunw_ggg,'(a)') window_str(:lnbc(window_str))
         close(lunw_ggg)
c
         if(dl.eq.char(92))then    ! Back-Slash (\)
            write(lunw_mul,'(a)') gggdir(:lrt)//'bin'//dl//'gfit '//
     $      filnam(:lnbc(filnam))
         else
            write(lunw_mul,'(a)') gggdir(:lrt)//'bin'//dl//'gfit '//
     $      filnam(:lnbc(filnam))//'>/dev/null'
         endif
      end do   !  kwin=1,mwin
      stop ' Increase parametr MWIN'
 103  close(lunr_win)
      close(lunw_mul)
c------------------------------------------------------------------
      nmod=0
      nvmr=0
c  Open the .mav and .ray files
      write(lunw_rpt,*)
c      write(86,*) 2,6
c      write(86,*) '  z    t    p    d    h2o  co2'
      open(lunw_mav,file=rlgfile(:lr-3)//'mav',status='unknown')
c  Open the file which will contain the slant paths
      write(lunw_mav,'(a)')version
      write(*,*)' Opening: ',rlgfile(:lr-3)//'ray'
      open(lunw_ray,file=rlgfile(:lr-3)//'ray',status='unknown')
      write(lunw_ray,'(i2,i4)') nlhead_ray,7+ncell+nlev
      write(lunw_ray,'(a)') version
      write(lunw_ray,'(a)') 'format='//ray_data_fmt
      write(ray_header_fmt,'(a3,i1,a9,i3.3,a10)')
     &  '(a,',ncell,'(a10,i1),',nlev,'(a8,i3.3))'
      write(lunw_ray,ray_header_fmt)
     &' SpectrumName                                                 '//
     &'   Zobs       Pobs       ASZA       Bend       FOV       Zmin  ',
     & ('     Cell_',j,j=1,ncell),
     & ('  Level_',j,j=0,nlev-1)
c      write(81,*) 2,7
c      write(81,*)'   obalt      asza      zenoff      solzen      urth  
c     &    zmin      bend'
c  Get the header/runlog information pertaining to RUNLAB
      write(6,*) ' Runlog=',gggdir(:lrt)//'runlogs'//dl//ext//dl//
     & rlgfile(:lr)
      open(lunr_rlg,file=gggdir(:lrt)//'runlogs'//dl//ext//dl//
     & rlgfile,status='old')
      call read_runlog_header(lunr_rlg,data_fmt_read_rl,col_labels_rl)

      do jspe=0,9999999
c1     continue
c      write(*,*)'Calling read_runlog_data...'
         call read_runlog_data_record(lunr_rlg,data_fmt_read_rl,
     &    col1,specname,iyr,idoy,zpdtim,oblat,
     &    oblon,obalt,asza,zenoff,azim,osds,
     &    opd,fovi,fovo,amal,ifirst,ilast,
     &    graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     &    tout,pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)
c         write(*,*)'Called read_runlog_data.'

c          write(*,*)'read_runlog_data_record:',istat,' '//specname(:20)
c         if(istat.ne.0) go to 99   !  EOF
c         if(col1.eq.':') go to 1
         if(istat.eq.3) stop 'Error reading runlog (e.g. "****")'
         if(istat.ne.0) exit   !  EOF
         if(col1.eq.':') cycle

cc  For satellite measurements, we want to use the tangent point
cc  location, not the observer location.
c         if(index(rlgfile,'atmos_merged').gt.0) then
c            call zenaz(2,oblat,oblon,obalt,iyr,1,idoy,zpdtim/24.0d0,
c     &      asza1,azim1,eorv,ervc,tplat,tplon,tpalt)
c            if((tplon-oblon).gt. 180) tplon=tplon-360.
c            if((tplon-oblon).lt.-180) tplon=tplon+360.
c            if(tplon.gt.360) tplon=tplon-360
c            if(tplon.lt.  0) tplon=tplon+360
c            oblat=tplat
c            oblon=tplon
c         endif

c         if(oblon.gt.180) oblon=oblon-360
c         if(oblon.lt.-180) oblon=oblon+360

         call julian(iyr,1,idoy,jd)
         trlg=dble(jd)+zpdtim/24

         call read_csi(specname,trlg,lunr_csi,ncell,tins,
     &   kcell,igas_in_cell,p_cell,t_cell,vmr_cell,cell_length)

c Check that specname(12:12) & cell_length are consistent
         if(modtype.eq.'FPIT') then
            if(specname(12:12).ne.'0' .and. cell_length(1).le.0.0)
     &      write(*,*)'Warning Inconsistent specname(12:12) & cell_info'
            if(specname(12:12).eq.'0' .and. cell_length(1).gt.0.0)
     &      write(*,*)'Warning Inconsistent specname(12:12) & cell_info'
         endif

         if(oblat.ge.0.0) then
            ns='N'
         else
            ns='S'
         endif
         if(oblon.ge.0.0) then
            ew='E'
         else
            ew='W'
         endif

         call julian(iyr,1,idoy,jul)
         frdoy=jul+zpdtim/24+oblon/360
         lspn=lnbc(specname)
         if(modtype.eq.'NCEP') then
c           Great Circle Distance (degrees)
            gcd=dacos(sin(d2r*oblat)*sin(d2r*latwas)+
     &      cos(d2r*oblat)*cos(d2r*latwas)*cos(d2r*(lonwas-oblon)))/d2r

c            diff=abs(int(frdoy)-int(frwas))+0.075*dabs(oblat-latwas)+
c     &       0.025*dabs(oblon-lonwas)
            if(ext.eq.'orb') then
               diff=4.8*abs(zpdtim-zpdwas)/24 + 0.042*gcd
            else
               diff=abs(int(frdoy)-int(frwas))+0.05*gcd
            endif
c            write(*,*) 'specname,frdoy,dd,gcd,diff = ', specname(:lspn),
c     &      frdoy,frdoy-frwas,gcd,diff

            if( diff.gt.0.6667) then
               call caldat(jul,iyyyy,imm,idd)
             write(modname_new,'(a5,i4.4,2i2.2,a1,i2.2,2a1,i3.3,a1,a4)')
     &        'NCEP_',iyyyy,imm,idd,'_',nint(abs(oblat)),ns,
     &         '_',nint(abs(oblon)),ew,'.mod'
               if(index(rlgfile,'atmos_merged').gt.0 .or.
     &            index(rlgfile,'at1ss').gt.0 .or.
     &            index(rlgfile,'at1sr').gt.0 .or.
     &            index(rlgfile,'at2ss').gt.0 .or.
     &            index(rlgfile,'at2ss').gt.0 .or.
     &            index(rlgfile,'at3ss').gt.0 .or.
     &            index(rlgfile,'at3sr').gt.0) then
                  modname_new=specname(:lspn)//'.mod'
c                  write(*,*) 'new modname = '//modname_new//specname
               endif
               zpdwas=zpdtim
               frwas=frdoy
               latwas=oblat
               lonwas=oblon
               jwas=jul
               vmrname_new=specname(1:lspn)//'.vmr' !  CIT/TCCON
c            vmrname_new=vmrname_bak ! ATMOS
            endif

         elseif(modtype.eq.'FPIT') then
            jhh = nhr*nint(zpdtim/nhr)
            ihh = mod(jhh+24,24)
            call caldat(jul+(jhh-ihh)/24,iyyyy,imm,idd)
            write(modname_new,
     &     '(a5,i4.4,3i2.2,a2,i2.2,2a1,i3.3,a1,a4)')
     &      'FPIT_',iyyyy,imm,idd,ihh,'Z_',nint(abs(oblat)),ns,
     &      '_',nint(abs(oblon)),ew,'.mod'
            write(vmrname_new,
     &     '(a4,i4.4,3i2.2,a2,i2.2,a1,a1,i3.3,a1,a4)')
     &      'JL1_',iyyyy,imm,idd,ihh,'Z_',nint(abs(oblat)),ns,
     &      '_',nint(abs(oblon)),ew,'.vmr'
         else
            stop ' unknown model type'
         endif
c         write(*,*) 'modname = '//modname_new
c         write(*,*) 'vmrname = '//vmrname_new
         frqcen=sngl((ifirst+ilast)*graw/2)
c
c  Fudge to allow nadir viewing (double path)
         if(asza .ge. 200.d0) then
            multipath=int(dabs(asza))/100
            asza=asza-100.0d0*multipath
         else
            multipath=1
         endif
         solzen=sngl(asza+zenoff)
         if(mod(jspe+1,1000).eq.0) write(*,*) jspe+1,' spectra....'

c  Following 25 lines removed 2020-06-07
c  Figure out whether the spectrum was through an internal cell.
c         do jcell=1,ncell
c            t_cell(jcell)=sngl(tins)+273.15
c            p_cell(jcell)=0.0
c            cell_length(jcell)=0.0
c            igas_in_cell(jcell)=0  
c         end do
c         if(lspn.ge.20) then
c            if(specname(1:4).ne.'InSb') then
c               if(specname(1:2).ne.'ss'
c     &         .and. specname(1:2).ne.'sr' 
c     &         .and. specname(1:3).ne.'ace' 
c     &         .and. specname(1:3).ne.'FTS') then ! skip ACE spectra
c                  ccell=specname(12:12)
c                  if( ccell.ne.'0') then
c                     if(ncell.ge.1) then
c                        cell_length(1)=0.0001 ! km (10cm)
c                        p_cell(1)=5.097       ! mbar ! TCCON Bruker internal cell
c                        iigas_in_cell(1)=15     ! HCl
c                        vmr_cell(1)=1.0       ! mole fraction
c                     endif  ! ncell.ge.1
c                  endif  ! ccell.ne.'0'
c               endif  ! specname(1:2)....
c            endif  ! specname(1:4)....
c         endif  ! lspn.ge.20
c
         if(specname(2:3).eq.'hg' .or. specname(2:3).eq.'in') then
c         if(specname(2:3).eq.'xg' .or. specname(2:3).eq.'xn') then
c            if(lspn.lt.8) then
c               write(*,*) 'lspn,specname=',lspn,specname
c               stop 'lspn<8'
c            endif
            if(ext.eq.'bal') then
               modname_new=ext//specname(4:8)//'.mod'  ! MkIV
               vmrname_new=ext//specname(4:8)//'.vmr'  ! MkIV
            elseif(ext.eq.'air')then
               modname_new='mod'//specname(4:)  ! MkIV
               vmrname_new='vmr'//specname(4:)  ! MkIV
            else                          ! ground-based MkIV
               if(modtype.eq.'NCEP') then
                  read(specname(4:8),'(i2,i3)')iyyyy,doy
c  Convert YYDDD-style date to YYYYMMDD
                  if(iyyyy.lt.80) then
                     iyyyy=iyyyy+2000
                  else
                     iyyyy=iyyyy+1900
                  endif
                  call julian(iyyyy,1,doy,jul)
                  call caldat(jul,iyyyy,imm,idd)
                  write(newmkivname,'(a2,i4,i2.2,i2.2)')
     &              rlgfile(1:2),iyyyy,imm,idd
                  vmrname_new=newmkivname//'.vmr'  ! MkIV
               endif
c
cc  Account for leaking OCS cell filled on day 279 of 2009
c               if(iyyyy.eq.2009) then
c                  igas_in_cell(1)=19       ! ocs
c                  pfill=1.37  ! Fill Pressure (mbar)
c                  p_cell(1)=pfill+0.060*(idoy-279+17)   ! Total Pressure (mbar).
c                  vmr_cell(1)=pfill/p_cell(1)
c               endif
c               if(  specname(4:16).eq.'09289.025_032' ! 20091016183900
c     &         .or. specname(4:16).eq.'09290.034_041' ! 20091017185600
c     &         .or. specname(4:16).eq.'09292.066_073' ! 20091019191800
c     &         .or. specname(4:16).eq.'09292.100_107' ! 20091019212000
c     &         .or. specname(4:16).eq.'09293.068_075' ! 20091020192000
c     &         .or. specname(4:16).eq.'09295.072_079' ! 20091022
c     &         .or. specname(4:16).eq.'09297.068_075' ! 20091024
c     &         .or. specname(4:16).eq.'09298.064_071' ! 20091025
c     &         .or. specname(4:16).eq.'09301.046_049' ! 20091028
c     &         .or. specname(4:16).eq.'09301.054_057' ! 20091028
c     &         .or. specname(4:16).eq.'09343.001_005' ! 20091209
c     &         .or. specname(4:16).eq.'09343.006_009' ! 20091209
c     &       ) then
c                  write(*,*) 'here ',specname(4:16)
c                  cell_length(1)=0.0001  ! km (=10 cm)
c               else
c                  cell_length(1)=0.00000  ! km 
c               endif
            endif

cc  NDSC HBr cell (no leaks)
c            if(ncell.ge.2) then
c               igas_in_cell(2)=16       ! HBr
c               p_cell(2)=2.25          ! Total Pressure (mbar)
c               vmr_cell(2)=0.8         ! mole fraction
c               if(  specname(4:16).eq.'09289.062_069'
c     &         .or. specname(4:16).eq.'09290.050_057'
c     &         .or. specname(4:16).eq.'09292.082_091'
c     &         .or. specname(4:16).eq.'09293.084_091'
c     &         .or. specname(4:16).eq.'09295.088_095'
c     &         .or. specname(4:16).eq.'09297.084_091'
c     &         .or. specname(4:16).eq.'09298.081_087'
c     &         .or. specname(4:16).eq.'09301.066_069'
c     &         .or. specname(4:16).eq.'15015.001_007'
c     &         .or. specname(4:16).eq.'17219.001_007'
c     &         .or. specname(4:16).eq.'17248.001_007'
c     &         .or. specname(4:16).eq.'19205.002_005'
c     &         .or. specname(4:16).eq.'19205.006_009'
c     &         .or. specname(4:16).eq.'19205.010_013'
c     &         .or. specname(4:16).eq.'19205.014_017'
c     &         .or. specname(4:16).eq.'19205.018_021'
c     &         .or. specname(4:16).eq.'19205.022_025'
c     &         .or. specname(4:16).eq.'19220.002_005'
c     &         .or. specname(4:16).eq.'19220.006_009'
c     &         .or. specname(4:16).eq.'19220.010_013'
c     &         .or. specname(4:16).eq.'19220.014_017'
c     &         .or. specname(4:16).eq.'19220.018_021'
c     &         .or. specname(4:16).eq.'19220.022_025'
c     &         .or. specname(4:16).eq.'19221.001_005'
c     &         .or. specname(4:16).eq.'19221.006_009'
c     &         .or. specname(4:16).eq.'19221.010_013'
c     &         .or. specname(4:16).eq.'19221.014_017'
c     &         .or. specname(4:16).eq.'19228.020_023'
c     &         .or. specname(4:16).eq.'19244.016_019'
c     &       ) then
c                  cell_length(2)=.00002  ! km (=2 cm)
c               else
c                  cell_length(2)=0.00000  ! km (2 cm)
c               endif
c            endif

c         elseif(specname(7:9).eq.'R0.') then  ! Kitt Peak
c            if(specname(1:1).eq.'7' .or. specname(1:1).eq.'8' .or.
c     &      specname(1:1).eq.'9' ) then
c            if(lspn.lt.8) stop 'lspn<8'
c              modname_new='kp19'//specname(1:6)//'.mod'
c              vmrname_new='kp19'//specname(1:6)//'.vmr'
c            else
c              modname_new='kp20'//specname(1:6)//'.mod'
c              vmrname_new='kp20'//specname(1:6)//'.vmr'
c            endif
         elseif(specname(1:2).eq.'ss' .or. specname(1:2).eq.'sr') then
            if(lspn.lt.6) stop 'lspn<6'
            sl = index(specname,'/')
            if(sl.eq.0)sl=7
            modname_new=specname(1:sl-1)//'.mod'    ! ACE
            vmrname_new=specname(1:sl-1)//'.vmr'    ! ACE
         elseif(specname(1:2).eq.'s1' .or. specname(1:2).eq.'r1') then
            if(lspn.lt.6) stop 'lspn<6'
            sl = index(specname,'/')
            if(sl.eq.0)sl=7
            modname_new=specname(1:sl-1)//'.mod'    ! ACE
            vmrname_new=specname(1:sl-1)//'.vmr'    ! ACE
         elseif(specname(1:24).eq.'FTS-TRANSMITTANCE-NA-NA-') then
            modname_new=rlgfile(1:lr-4)//'.mod'    ! ACE (new names)
            vmrname_new=rlgfile(1:lr-4)//'.vmr'    ! ACE (new names)
         elseif(rlgfile(lr-1:lr).eq.'ws') then
            if(lspn.lt.6) stop 'lspn<6'
            modname_new=specname(2:6)//'.zpt' ! Wollongong spectra
            vmrname_new=specname(2:6)//'.ref' ! Wollongong spectra
c         elseif(specname(1:2).eq.'iz') then
c            if(lspn.lt.6) stop 'lspn<6'
c            modname_new='pt_'//specname(1:10)//'.prf' ! Izana
c            vmrname_new='vmr_'//specname(1:10)//'.ref' ! Izana
         elseif(specname(1:2).eq.'o2') then
            if(lspn.lt.6) stop 'lspn<6'
            modname_new=specname(1:6)//'av.mod' ! TMF FTUVS O2 A-band spectra
            vmrname_new=specname(1:6)//'av.vmr' ! TMF FTUVS O2 A-band spectra
c         elseif(specname(1:1).eq.'9' .or. specname(1:1).eq.'0') then
c            if(lspn.lt.6) stop 'lspn<6'
c            modname_new=specname(1:6)//'.mod'  ! Kitt Peak
c            vmrname_new=specname(1:6)//'.vmr'  ! Kitt Peak
         elseif(specname(1:3).eq.'j60' ) then
            modname_new=specname(1:lspn)//'.mod'
            vmrname_new=specname(1:lspn)//'.vmr'
c         elseif(ext.eq.'orb')then
c            modname_new='mal'//specname(4:10)//'.mod'
c            vmrname_new='mal'//specname(4:10)//'.vmr'
         else
c            write(*,*) lspn,specname(:lspn)
c            if(lspn.lt.10) stop 'lspn<10'
c            modname_new=specname(1:10)//'.mod' !  CIT/TCCON
c            vmrname_new=specname(1:10)//'.vmr' !  CIT/TCCON
         endif
c         lm=lnbc(modname_new)
c         lv=lnbc(vmrname_new)

c         write(*,*)' modname_new =  ', modname_new(:32)
c         write(*,*)' modname_cur =  ', modname_cur(:32)
         newmod=.false.  ! Current model okay
c         if(modname_new(:lm).ne.modname_cur(:lm)) then 
         if(modname_new.ne.modname_cur) then 
            inquire(file=gggdir(:lrt)//'models'//dl//ext//dl//
     &      modname_new,exist=newmod_exist)
            ninqmod=ninqmod+1
            if(newmod_exist .eqv. .true.) then
               modname_cur=modname_new
               nmod=nmod+1
            else
               write(6,'(a)')'specname = '//specname(:50)
               write(6,'(a)')'Didnt find model: '//modname_new(:40)
               write(lunw_rpt,'(a)')'Didnt find model:',modname_new(:50)
               if(modtype.eq.'FPIT') stop 'Didnt find model'
               modname_cur=modname_bak
            endif
            newmod=.true.
            date_mod=iyr+idoy/365.25 !from the runlog
            write(*,'(a)')'Calling read_model_fc: '//
     &      gggdir(:lrt)//'models'//dl//ext//dl//modname_cur(:32)
            call read_model_fc(lunr_mod,gggdir(:lrt)//'models'//dl//
     &      ext//dl//modname_cur,z,mmw,nlev,t,p,d,h2o_dmf,roc,tlat)

c Compute roc in the E-W direction (since we dont have the AZIM)
            rocx=roc
c            if(index(rlgfile,'atmos_merged.orl').gt.0 .or.
c     &      index(rlgfile,'at1ss07.orl').gt.0 ) then
c            rocx=roc/sqrt(cos(d2r*oblat)**2+(0.99665*sin(d2r*oblat))**2)
c            endif

            if(nlev.gt.1) then
               ztrop_gct=compute_ztrop_gct(nlev,z,t)
               ztrop_mod=compute_ztrop(nlev,z,t)
            else
               ztrop_mod=0.0
            endif
         endif  !  (modname_new.ne.modname_cur)

c   Select and read .vmr file
c         write(*,*)' vmrname_new =  ', vmrname_new(:20)
c         write(*,*)' vmrname_cur =  ', vmrname_cur(:20)
         newvmr=.false.
c         if(vmrname_new(:lv).ne.vmrname_cur(:lv)) then ! search for new vmr file
         if(vmrname_new.ne.vmrname_cur) then ! search for new vmr file
            inquire(file=gggdir(:lrt)//'vmrs'//dl//ext//dl//vmrname_new,
     &      exist=newvmr_exist)
            ninqvmr=ninqvmr+1
            if(newvmr_exist .eqv. .true.) then
               vmrname_cur=vmrname_new
               nvmr=nvmr+1
               newvmr=.true.
            else
c               lv=lnbc(vmrname_new)
               write(lunw_rpt,*)'Didnt find vmrfile: '//vmrname_new(:20)
     &         //modtype
               if(modtype.eq.'FPIT') stop 'Missing vmr file.'
               vmrname_cur=vmrname_bak
               write(lunw_rpt,*)'So instead using: '//vmrname_cur(:20)
            endif

c            newvmr=.true.
            write(*,'(a)')'Calling read_refvmrs: '//
     &      gggdir(:lrt)//'vmrs'//dl//ext//dl//vmrname_cur(:28)

            call read_refvmrs(lunr_vmr,
     &      gggdir(:lrt)//'vmrs'//dl//ext//dl//vmrname_cur,nlev,
     &      z,mgas,gggdir(:lrt)//'models'//dl//ext//dl//modname_new,
     &      vmrlabel,refvmr,ngas,strend,gradlat,seacycle,reflat_vmr,
     &      date_vmr,ztrop_vmr)

c If there's only 1 level (i.e. lab), or the desired vmr file for that
c particular measurement already exists, don't adjust the .vmr file
            if(nlev.gt.1.and.modtype.ne.'FPIT'.and. ext.ne.'bal') then
               write(*,*)'Adjusting VMRs',date_mod,date_vmr
               call calc_itcz(oblon,imm,itcz_lat,itcz_width)
               call resample_vmrs_at_effective_altitudes(nlev,z,mgas,
     &         ngas,itcz_lat,itcz_width,
     &         refvmr,ztrop_mod,ztrop_vmr,oblat,reflat_vmr,apvmr)
               call apply_vmr_latitude_gradients(nlev,z,mgas,ngas,
     &         gradlat,apvmr,ztrop_mod,reflat_vmr,oblat,apvmr)
               call apply_secular_trends(nlev,z,mgas,ngas,strend,
     &         apvmr,ztrop_mod,reflat_vmr,oblat,date_mod,date_vmr,
     &         apvmr)
               fryr=date_mod-int(date_mod)
               call apply_seasonal_cycle(nlev,z,mgas,ngas,seacycle,
     &         ztrop_mod,oblat,fryr,apvmr)
c               write(*,*) 'co2 apvmr = ',(apvmr(2,j),j=1,9)
            else
               write(*,*)'Setting VMR to refvmr (no adjustment)'
               do ilev=1,nlev
                  do jgas=1,mgas
                     apvmr(jgas,ilev)=refvmr(jgas,ilev)
                  end do
               end do
            endif
         endif   !  if(vmrname_new.eq.vmrname_cur) then

c  Output model information (SUNRUN.MAV)
         if(newvmr.or.newmod) then
            write(*,*)'Writing to .mav: Next Spectrum: '//specname
            write(lunw_rpt,*)'Writing to mav: Next Spectrum: '//specname
            t_cell(kcell)=sngl(tins)+273.15
            write(lunw_mav,'(a)')'Next Spectrum:'//specname
            write(lunw_mav,'(3i4)') 6, nspeci+4, nlev+ncell
            write(lunw_mav,'(a,f9.3)') 'Tropopause Altitude:',ztrop_mod
            write(lunw_mav,'(a,f9.3)') 'Observer Latitude:',oblat
            write(lunw_mav,'(a)') gggdir(:lrt)//'vmrs'//dl//ext//dl//
     &      vmrname_cur(:lnbc(vmrname_cur))
            write(lunw_mav,'(a)') gggdir(:lrt)//'models'//dl//ext//dl//
     &      modname_cur(:lnbc(modname_cur))
            zpbl=0.0
c            write(*,*)'newvmr_exist: ',newvmr_exist,vmrname_new(:lv)
c            if (newvmr_exist) then
c              Do nothing when the vmr file is found.
c            else
c  For ground-based observations, over-write apriori H2O vmr with model.
c  For HDO, the factor 0.14*(8.0+log10(h2o_dmf(ilev))) adjusts for the
c  effects of isotopic fractionation (Rayleigh Distillation curve).
c  For H2O = 0.03,  HDO=0.91*H2O (sea-level in tropics)
c  For H2O = 0.01,  HDO=0.84*H2O (sea-level at mid-latitude)
c  For H2O = 3E-06, HDO=0.35*H2O (tropopause)
c  Changed from 0.16 on 2019-11-27 on the basis of TCCON data
               if(abs(h2o_dmf(1)).gt.0.0) then ! An H2O profile was found in the .mod file
                  if ( pout.gt.600.0 ) then ! ground-based observations
                     write(*,*)'Replacing H2O & HDO vmrs with '//modtype
                     do ilev=1,nlev
                        apvmr(1,ilev)=abs(h2o_dmf(ilev))  ! H2O
                        apvmr(49,ilev)=abs(h2o_dmf(ilev))*
     &                  0.14*(8.0+log10(abs(h2o_dmf(ilev)))) ! HDO
                     end do
                  endif
               endif
c            endif  !  if (newvmr_exist) then

c            write(*,*) 'write_mav: t_cell',(t_cell(j),j=1,ncell)
c            write(*,*) 'write_mav: p_cell',(p_cell(j),j=1,ncell)
c            write(*,*) 'write_mav: vmr_cell',(vmr_cell(j),j=1,ncell)
c            write(*,*) 'write_mav: igas_in_cell',(igas_in_cell(j),j=1,ncell)
c      Level 3 corresponds to 2 km, the average altitude of MLO at SMO
c      Gases 2,4,6 are CO2, N2O, CH4
c            write(*,*)'Calling write_mav....'
            call write_mav(z,t,p,d,apvmr,nlev,lunw_mav,ncell,
     &      t_cell,p_cell,vmr_cell,igas_in_cell,vmrlabel,isofile)
c            write(*,*)'Called write_mav.'
         endif
c------------------------------------------------------------------

         if(pout.le.0.0) then
            zpres=sngl(obalt)
         else
            zpres=height(sngl(pout)/1013.25,p,t,z,nlev)
         endif
c
         if(abs(sngl(obalt)-zpres).gt.0.1) write(lunw_rpt,*)
     $   'Warning: Geometrical & Pressure altitudes differ:',obalt,zpres
c
c  Determine ALT from POBS if > 400.0    Otherwise use geometric altitude
         if(pout.lt.400.0d0 .and. iyr.lt.1979) then
            zobs=sngl(obalt)
            stop ' pout.lt.400.0d0 .and. iyr.lt.1979'
         else
            zobs=zpres   ! prefer to determine GRADIUS from POBS
         endif
c
c  Determine ray slant paths
         fovr=90.*sngl(fovo)/spi     ! convert radians diameter to deg radius
c         write(*,*)'Calling TLPATH...',rocx,zobs,solzen,wavtkr,frqcen
         call tlpath(nlev,z,t,p,solzen,fovr,rocx,zobs,
     &    sngl(wavtkr),frqcen,zmin,bend,splos,rc)
c         write(*,*)'Called TLPATH: zmin,bend =',zmin,bend
c         write(81,'(7f11.3)') obalt,asza,zenoff,solzen,
c     &    (6378.+obalt)*sin(d2r*solzen)-6378.0,zmin,bend
         if(rc.ne.0) write(6,'(a32,5f10.4)')'Error: TLPATH:'//specname,
     &   zobs,asza,zenoff,solzen,zmin
c
c  Distribute AIPL between the two atmospheric levels that bracket POUT
         do ilev=2,nlev
            if(pout/1013.25.gt.p(ilev)) exit
         end do
         if(ilev.gt.nlev) then  ! Failed to find low enough pressure level
            if(aipl.gt.0) then
               write(*,*) 'Warning: aipl>0; pout<p(i)'
               write(*,*) specname(:20), aipl,pout,p(nlev),nlev
            endif
c            if(nlev.gt.1) stop 'vmrs dont extend high enough'
            splos(nlev)=splos(nlev)+sngl(aipl*(pout/1013.25))/p(nlev)
            if(nlev.eq.1) splos(nlev)=sngl(aipl)
         else
            delta_p=p(ilev-1)-p(ilev)
            splos(ilev-1)=splos(ilev-1)+
     &                    sngl(aipl*(pout/1013.25-p(ilev)))/delta_p
            splos(ilev)=splos(ilev)+
     &                  sngl(aipl*(p(ilev-1)-pout/1013.25))/delta_p
         endif
c
         write(lunw_rpt,'(a)') 'Writing to .ray: '//specname(:40)
c         write(96,*) specname(:20),kcell,cell_length(1),cell_length(2)
         write(lunw_ray,ray_data_fmt)specname,obalt,pout,
     &    solzen,bend,fovo,zmin,(cell_length(j),j=1,ncell),
     &    (wlimit(dble(multipath)*splos(j),'f10.5'),j=1,nlev)
      end do  ! jspe = 0,999999

c      go to 1
c 99   close(lunr_rlg)
      close(lunr_rlg)
      close(lunw_ray)
      close(lunw_mav)
      close(lunw_rpt)
      if(jspe.le.0)write(6,*)'none of these spectra could be accessed'
      if(jspe.le.0)write(6,*)'or error reading runlog (wrong format?)'
      write(*,*) jspe,' valid spectra found in runlog'
      write(*,*) nmod,' new .mod files read from',ninqmod,' inquires'
      write(*,*) nvmr,' new .vmr files read from',ninqvmr,' inquires'

c  Code to generate post_processing.sh batch file and associated inputs.

      call write_postprocessfile(ext,rlgfile,modtype)

c      do ih=1,100
c        write(84,*) ih, dhp*ih, histop(ih)
c        write(85,*) ih, dhz*ih, histoz(ih)
c      end do

      stop
      end
