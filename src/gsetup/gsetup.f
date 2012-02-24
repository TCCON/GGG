c   gsetup.f
c   Provides user interface to create input files for gfit
c
c  Inputs:
c     Menu-driven user input
c
c  Outputs:
c     Multigg.bat  !  Batch file
c     gas_1234_runlog.ggg  ! GFIT input files (one for each window)
c     runlog.mav           ! P/T/vmr profiles interpolated onto chosen levels
c     runlog.ray           ! Ray-traced atmospheric slant paths
c     gsetup.rpt           ! File containing report

c  If GSETUP asks you for a model or vmr, it couldn't find a
c  model or vmr in the appropriate sub-directory having the default name.
c

      implicit none
      include "../ggg_const_params.f"
      include "../ggg_int_params.f"
      include "int_params.f"

      integer*4 apo,bytepw,fbc,fnbc,ifirst,ilast,ilev,jcell,
     & interp,idoy,iyr,j,k,w1,w2,le,lg,lc,lnbc,lrt,
     & lunw_ray,lunw_ggg,
     & lunr_rlg,lunw_mul,
     & lunw_mav,lunw_rpt,lunr_mod,
     & lunr_vmr,
     & lunr_iso,lunw_iso,
     & iday,jgas,
     & nlev,nlhead_rlg,nlhead_ray,nlhead_ggg,nspeci,ngas,
     & lr,lspn,i,ncol,doy,iyyyy,imm,idd,jul,
     & multipath, nspe,possp,rc,istat,lm,lv,
     & lmn,lvn,lunr_men
      parameter (lunw_ggg=56)    ! for writing the .ggg files
      parameter (lunr_iso=57)    ! for reading isotopologs.dat
      parameter (lunw_iso=58)    ! for writing isotopologs_local.dat
      parameter (lunw_mav=60)    ! for writing the .mav file
      parameter (lunr_men=61)    ! for reading the various .men files
      parameter (lunr_mod=63)    ! for reading the .mod file (readmodFC)
      parameter (lunw_mul=66)    ! for writing the multiggg.sh file
      parameter (lunw_ray=68)    ! for writing the .ray file
      parameter (lunr_rlg=69)    ! for reading the runlog
      parameter (lunw_rpt=72)    ! for writing gsetup.rpt
      parameter (lunr_vmr=73)    ! for reading the .vmr file (readvmrFC)
      parameter (nlhead_ray=3)   ! Number of header lines in .ray file.
      parameter (nlhead_ggg=18)  ! Number of header lines in .ggg file.
      real*4 erroff
      parameter (erroff=0.004,apo=2,interp=1)
      logical*4 newmod,newvmr,vmr_found

      integer*4
     & gas_in_cell(ncell)

      real*4 dum,freq,frqcen,height,tlat,
     & zmin,zpres,compute_seasonal_cycle,
     & p_cell(ncell),t_cell(ncell),vmr_cell(ncell),
     & zobs,             ! the geocentric radius of the observer fed to TLPATH
     & bend,             ! total ray bending due to refraction
     & d(mlev),          ! number densities at levels (molec.cm-3)
     & h2ovmr(mlev),     ! h2o profile from model file
c     & co2vmr(mlev),     ! co2 profile from simulate_co2_vmr
     & fovr,             ! External angular radius of FOV in degrees
     & mmw(mlev),        ! mean molecular weights
     & p(mlev),          ! pressures of levels (atm.)
     & r(mlev),          ! radii of levels (km) = z + roc
     & roc,              ! radius of curvature (km)
     & cell_length(ncell),      ! Length of cell permanently in solar beam
     & solzen,           ! solar zenith angle used (may be refracted)
     & splos(mlev),      ! array of line-of-sight slant paths
     & t(mlev),          ! temperatures of levels (K)
     & refvmr(mgas,mlev), ! buffer for reference vmr's
     & apvmr(mgas,mlev), ! buffer for a priori vmr's
     & z(mlev)           ! altitudes of levels (km)

      real*8  fryr,z8,
     & reflat_vmr,date_mod,date_vmr,ztrop_vmr

      real*8  zpbl, ztrop_ncep,ztrop_gct, delta_p, wlimit,
     & oblat,            ! observation latitude (deg).
     & oblon,            ! observation longitude (deg).
     & obalt,            ! observation altitude (km)
     & asza,             ! astronomical solar zenith angle (unrefracted)
     & opd,              ! Optical path difference (cm) of interferogram
     & graw,             ! spacing of raw spectrum (cm-1) from GETINFO
     & zpdtim,           ! Time of ZPD (UT hours)
     & zenoff,           ! Zenith angle pointing offset (deg)
     & azim,             ! Solar Azimuth Angle
     & osds,             ! Observer-Sun Doppler Stretch (ppm_
     & fovi,             ! Internal angular diameter of FOV (radians)
     & fovo,             ! External angular diameter of FOV (radians)
     & amal,             ! angular misalignment of interferometer (radians)
     & zoff,             ! Zero level offset (dimensionless fraction)
     & snr,              ! Signal-to-Noise Ratio (dimensionless)
     & tins,             ! Inside temperature
     & pins,             ! Inside pressure
     & hins,             ! Inside humidity
     & tout,             ! Outside temperature
     & pout,             ! Outside pressure
     & hout,             ! Outside humidity
     & sia,              ! Solar Intensity Average (arbitrary units)
     & fvsi,             ! Fractional Variation in Solar Intensity
     & wspd,             ! Wind Speed (m/s)
     & wdir,             ! Wind Direction (deg)
     & aipl,             ! Airmass-Independent Path length (km)
     & lasf,             ! Laser Frequency (e.g. 15798 cm-1)
     & wavtkr            ! Frequency of the sun-tracker (actively-tracked)

        
      character
     & window*140,         ! name of window
     & filnam*48,          ! general purpose character buffer
     & filnamwas*48,       ! general purpose character buffer
     & header*92,          ! general purpose character buffer
     & levels*24,          ! name of file containing levels
     & listof*24,          ! the selected list of windows
     & newmkivname*10,     ! path to new model
     & modname*80,         ! path to new model
     & newmodname*80,      ! path to new model
     & vmrname*80,         ! path to new vmr
     & newvmrname*80,      ! path to new vmr
     & menuinput*(mfilepath),       ! path to xxx.men file
     & prvwin*14,          ! name of previous window (i.e. GAS_1234)
     & slk*1,              ! symbol which links gas name and frequency
     & ccell*1,            ! Character in specname denoting internal cell 
c    & user*8,             ! investigator
     & ray_fmt*26,
     & isofile*(mfilepath),         ! path to "isotopolog.dat"
     & vmrlabel*1024,      ! column labels from vmr file
     & str_isotop*294      ! line from isotopologs.dat file
      character
     & col1*1,             !first column of runlog record
     & apf*2,              !apodization function (e.g. BX N2, etc)
     & dl*1,               !forward or backward slash
     & ext*3,              !geometry ['air','bal','gnd','lab',orb','syn']
     & gggdir*(mpath),     !ggg directory path (GGGPATH?)
     & specname*(nchar),   !spectrum name
     & version*64,         !current program version
     & rlgfile*120         !name of runlog file


      version=
     & ' GSETUP                   Version 3.4.2    05-Nov-2011    GCT '
      modname='                                                '
      vmrname='                                                '
      filnamwas='qwertyuioqwertyuioqwertyuioqwertyuio'
      vmr_found=.false.

c     Platform specification:      DG090519
      call get_ggg_environment(gggdir, dl)
      lrt=lnbc(gggdir)     !Length of gggdir
c
c  Find NSPECI and copy isotopologs.dat to the local directory
      isofile=gggdir(:lrt)//'isotopologs/isotopologs.dat'
      open(lunr_iso,file=isofile,status='old')
      open(lunw_iso,file='./isotopologs_local.dat',status='unknown')
      do i=1,999
         read(lunr_iso,'(a)',end=76) str_isotop
         write(lunw_iso,'(a)') str_isotop(:lnbc(str_isotop))
      end do
76    nspeci=i-1
      close(lunr_iso)
      close(lunw_iso)

      slk='_'
      ext(1:3)='   '
      open(lunw_rpt,file='gsetup.rpt', status='unknown')

c  choose an observation geometry
      write(6,*) version
 3    write(6,9913)
 9913 format(' Geometry (a=air,b=bal,g=gnd,l=lab,o=orb,s=syn) ? ',$)
      read(5,'(a)') ext(1:1)
      if(ext(1:1).eq.'a') ext(1:3)='air'
      if(ext(1:1).eq.'b') ext(1:3)='bal'
      if(ext(1:1).eq.'g') ext(1:3)='gnd'
      if(ext(1:1).eq.'l') ext(1:3)='lab'
      if(ext(1:1).eq.'o') ext(1:3)='orb'
      if(ext(1:1).eq.'s') ext(1:3)='syn'
      if(ext(2:3).eq.'  ') go to 3

c      if(ext(1:1).eq.'b' .or. ext(1:1).eq.'b') then
c          og='l'
c      else
c          og='v'
c      endif
c------------------------------------------------------------------
c  Choose which runlog to analyze
      menuinput=gggdir(:lrt)//'runlogs'//dl//ext//dl//'runlogs.men'
      call readmenu(lunr_men,menuinput,rlgfile)
      lr=lnbc(rlgfile)
c------------------------------------------------------------------
c  choose which levels to use
      menuinput=gggdir(:lrt)//'levels'//dl//'levels.men'
      call readmenu(lunr_men,menuinput,levels)
c  Read chosen level altitudes
      open(lunr,file=gggdir(:lrt)//'levels'//dl//levels(:lnbc(levels)),
     $status='old')
      do ilev=1,mlev
         read(lunr,*,end=877) z(ilev),mmw(ilev)
      end do
      read(lunr,*,end=877) dum
      write(6,*) 'Warning: number of levels exceeds NLMAX'
 877  nlev=ilev-1
      close(lunr)
      if(nlev+1.gt.mlev) stop 'Error: nlev+1 > mlev '
c------------------------------------------------------------------
c  choose which list of windows to analyze
      menuinput=gggdir(:lrt)//'windows'//dl//ext//dl//'windows.men'
      call readmenu(lunr_men,menuinput,listof)

c  read the individual microwindows and create the xxxxxxxx.ggg input files
      if(dl.eq.char(92))then  ! Back-Slash (\)
          open(lunw_mul,file='multiggg.bat',status='unknown')
      else
          open(lunw_mul,file='multiggg.sh',status='unknown')
      endif
      open(lunr,file=gggdir(:lrt)//'windows'//dl//ext
     $ //dl//listof(:lnbc(listof)),status='old')
      read(lunr,'(a91)') header
      prvwin='            '
      do k=1,999999
 666    read(lunr,'(a)',end=103) window
        if(window(1:1).eq.':' .or. window(1:1).eq.';') go to 666
        w1=fnbc(window)
        if(w1.eq.0) go to 666 ! blank line
        w2=index(window,'.')
c        write(*,*) w1,w2,window(w1:w2)
        read(window(w1:w2),*) freq
        write(6,'(a)') ' '//window(:lnbc(window))

        lc=index(window,':')  ! lc=position of colon
        lg=lc+fnbc(window(lc+1:))
        le=lg+fbc(window(lg+1:))-1
        filnam=window(lg:le)//slk//window(w1:w2)//rlgfile(:lr-3)//'ggg'
        if(filnam .eq. filnamwas) then
c          slk=char(ichar(slk)-1) ! change to the previous ASCII character
          slk='-'
        else
          slk='_'
        endif
        filnam(le-lg+2:le-lg+2)=slk
c        write(*,*)lg,le,filnam
        filnamwas=filnam
c        write(*,*)'filnam=',filnam
        open(lunw_ggg,file=filnam,status='unknown')
        write(lunw_ggg,'(i3)') nlhead_ggg
        write(lunw_ggg,'(a)') version
c        write(lunw_ggg,'(a)') version
        write(lunw_ggg,'(a)')
     &      gggdir(:lrt)//'config'//dl//'data_part.lst'
        write(lunw_ggg,'(a)')
     &      gggdir(:lrt)//'apriori'//dl//'gfit_ap.'//ext
        write(lunw_ggg,'(a)') gggdir(:lrt)//'runlogs'//dl//ext//dl//
     $  rlgfile(:lr)
        write(lunw_ggg,'(a)') gggdir(:lrt)//'levels'//dl//levels
        write(lunw_ggg,'(a)') gggdir(:lrt)//'models'//dl//ext//dl
        write(lunw_ggg,'(a)') gggdir(:lrt)//'vmrs'//dl//ext//dl
        write(lunw_ggg,'(a)') rlgfile(:lr-3)//'mav'
        write(lunw_ggg,'(a)') rlgfile(:lr-3)//'ray'
        write(lunw_ggg,'(a)') 'isotopologs_local.dat'
        write(lunw_ggg,'(a)')
     &      gggdir(:lrt)//'windows'//dl//ext//dl//listof

c  write names of linelists
c        write(lunw_ggg,'(a)')      ! llsize
        write(lunw_ggg,'(a)')
     &  gggdir(:lrt)//'linelist'//dl//'atm.101 '//
     &  gggdir(:lrt)//'linelist'//dl//'gct.101 '//
     &  gggdir(:lrt)//'linelist'//dl//'fcia.101 '//
     &  gggdir(:lrt)//'linelist'//dl//'scia.101'
        write(lunw_ggg,'(a)')
     &  gggdir(:lrt)//'linelist'//dl//'solar_merged.108'

c  write location of .ak files (averaging kernels)
        write(lunw_ggg,'(a)') gggdir(:lrt)//'ak'//dl//'k'

c  write location of spectral fits
        write(lunw_ggg,'(a)') gggdir(:lrt)//'spt'//dl//'z'

        write(lunw_ggg,'(a)') filnam(:lnbc(filnam)-3)//'col'
c        write(lunw_ggg,'(a)') filnam(lnbc(filnam)-11:lnbc(filnam)-3)//'col'
        write(lunw_ggg,'(a)') window
        close(lunw_ggg)
c
        if(dl.eq.char(92))then    ! Back-Slash (\)
          write(lunw_mul,'(a)') gggdir(:lrt)//'bin'//dl//'gfit<'//
     $    filnam(:lnbc(filnam))
        else
          write(lunw_mul,'(a)') gggdir(:lrt)//'bin'//dl//'gfit<'//
     $    filnam(:lnbc(filnam))//'>/dev/null'
        endif
      end do   !  k=1,999999
 103  close(lunr)
      close(lunw_mul)
c------------------------------------------------------------------
c  Compute the slant paths and write them to disk.
      write(lunw_rpt,*)
      write(lunw_rpt,'(a)')
     $'  Spectrum        Zobs   Pobs    ASZA    BEND    FOVO    TANG'
      write(lunw_rpt,'(a)')
     $'                   km    mbar    deg.    deg.    mrad     km'
      nspe=0
      open(lunw_mav,file=rlgfile(:lr-3)//'mav',status='unknown')
c  Open the file which will contain the slant paths
      write(lunw_mav,'(a)')version
      write(*,*)' Opening: ',rlgfile(:lr-3)//'ray'
      open(lunw_ray,file=rlgfile(:lr-3)//'ray',status='unknown')
      write(lunw_ray,'(i2,i4)') nlhead_ray,7+ncell+nlev
      write(lunw_ray,'(a)') version
c      ray_fmt='(a,2(a5,i1),250(a2,i3.3))'
      write(ray_fmt,'(a3,i1,a9,i3.3,a10)')
     &  '(a,',ncell,'(a10,i1),',nlev,'(a8,i3.3))'
      write(lunw_ray,ray_fmt)
     &' SpectrumName                           '//
     &'   Zobs       Pobs       ASZA       Bend       FOV       Zmin  ',
     & ('     Cell_',jcell,jcell=1,ncell),
     & ('  Level_',j,j=0,nlev-1)
c  Get the header/runlog information pertaining to RUNLAB
      write(6,*) ' Runlog=',gggdir(:lrt)//'runlogs'//dl//ext//dl//
     & rlgfile(:lr)
      open(lunr_rlg,file=gggdir(:lrt)//'runlogs'//dl//ext//dl//
     & rlgfile,status='old')
      read(lunr_rlg,*,err=888) nlhead_rlg,ncol           ! read runlog header line
      do i=2,nlhead_rlg
         read(lunr_rlg,*)
      end do
888   continue  !  supports older runlog formats
1     call read_runlog(lunr_rlg,col1,specname,iyr,idoy,zpdtim,oblat,
     & oblon,obalt,asza,zenoff,azim,osds,
     & opd,fovi,fovo,amal,ifirst,ilast,
     & graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     & tout,pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)
       lspn=lnbc(specname)
c       write(*,*)'read_runlog:',istat,la,' '//specname(:20)
       if(istat.ne.0) go to 99   !  EOF
       if(col1.eq.':') go to 1
       frqcen=(ifirst+ilast)*graw/2
c
c  Fudge to allow nadir viewing (double path)
      if(asza .ge. 200.d0) then
        multipath=int(dabs(asza))/100
        asza=asza-100.0d0*multipath
      else
        multipath=1
      endif
      solzen=sngl(asza+zenoff)
      nspe=nspe+1
      if(mod(nspe,1000).eq.0) write(*,*) nspe,' spectra....'

c  Figure out whether the spectrum was through an internal cell.
      do jcell=1,ncell
         t_cell(jcell)=tins+273.16
         p_cell(jcell)=0.0
         cell_length(jcell)=0.0
         gas_in_cell(jcell)=0  
      end do
      if(lspn.ge.20) then
         if(specname(1:4).ne.'InSb') then
         if(specname(1:2).ne.'ss'
     &   .and. specname(1:2).ne.'sr' 
     &   .and. specname(1:3).ne.'FTS') then ! skip ACE spectra
            ccell=specname(12:12)
            if( ccell.ne.'0') then
               if(ncell.ge.1) then
                  cell_length(1)=0.0001 ! km (10cm)
                  p_cell(1)=5.097       ! mbar ! TCCON Bruker internal cell
                  gas_in_cell(1)=15     ! HCl
                  vmr_cell(1)=1.0       ! mole fraction
               endif  ! ncell.ge.1
            endif  ! ccell.ne.'0'
         endif  ! specname(1:2)....
         endif  ! specname(1:4)....
      endif  ! lspn.ge.20
c
      if(specname(2:3).eq.'hg' .or. specname(2:3).eq.'in') then
         if(lspn.lt.8) stop 'lspn<8'
c         newmodname=rlgfile(1:5)//specname(4:8)//'.mod'  ! MkIV
c         newvmrname=rlgfile(1:5)//specname(4:8)//'.vmr'  ! MkIV
         if(ext.eq.'bal') then
            newmodname=ext//specname(4:8)//'.mod'  ! MkIV
            newvmrname=ext//specname(4:8)//'.vmr'  ! MkIV
         elseif(ext.eq.'air')then
            newmodname='mod'//specname(4:)  ! MkIV
            newvmrname='vmr'//specname(4:)  ! MkIV
c            cell_length(1)=0.002  ! km (1m)
c            p_cell(1)=655.0       ! mbar !
c            gas_in_cell(1)=1      ! H2O
c            vmr_cell(1)=-0.02     ! mole fraction
         else
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
     &        rlgfile(1:2),iyyyy,imm,idd
         newmodname=newmkivname//'.mod'  ! MkIV
         newvmrname=newmkivname//'.vmr'  ! MkIV
         gas_in_cell(1)=19       ! ocs
         cell_length(1)=0.00000  ! km (10 cm)
         p_cell(1)=1.37+0.08*(idoy-279)   ! Total Pressure (mbar). OCS cell filled on day 279
         vmr_cell(1)=1.37/p_cell(1)
         if(specname(4:16).eq.'09289.025_032'
     &   .or. specname(4:16).eq.'09292.066_073'
     &   .or. specname(4:16).eq.'09292.100_107'
     &   .or. specname(4:16).eq.'09293.068_075'
     &   .or. specname(4:16).eq.'09295.072_079'
     &   .or. specname(4:16).eq.'09297.068_075'
     &   .or. specname(4:16).eq.'09298.064_071'
     &   .or. specname(4:16).eq.'09301.046_049'
     &   .or. specname(4:16).eq.'09301.054_057'
     &   .or. specname(4:16).eq.'09343.001_005'
     &   .or. specname(4:16).eq.'09343.006.009') cell_length(1)=0.00010
         endif
         if(ncell.ge.2) then
         gas_in_cell(2)=16       ! HBr
         cell_length(2)=0.00000  ! km (2 cm)
         p_cell(2)=2.25          ! Total Pressure (mbar)
         vmr_cell(2)=0.8         ! mole fraction
         if(specname(4:16).eq.'09289.062_069'
     &   .or. specname(4:16).eq.'09292.082_091'
     &   .or. specname(4:16).eq.'09293.084_091'
     &   .or. specname(4:16).eq.'09295.088_095'
     &   .or. specname(4:16).eq.'09297.084_091'
     &   .or. specname(4:16).eq.'09298.081_087'
     &   .or. specname(4:16).eq.'09301.066_069'
     &   .or. specname(4:16).eq.'09301.066_069') cell_length(2)=0.00002
         endif
      elseif(specname(7:9).eq.'R0.') then  ! Kitt Peak
         if(specname(1:1).eq.'7' .or. specname(1:1).eq.'8' .or.
     &   specname(1:1).eq.'9' ) then
         if(lspn.lt.8) stop 'lspn<8'
c           newmodname=specname(1:8)//'.mod'
c           newvmrname=specname(1:8)//'.vmr'
           newmodname='kp19'//specname(1:6)//'.mod'
           newvmrname='kp19'//specname(1:6)//'.vmr'
         else
           newmodname='kp20'//specname(1:6)//'.mod'
           newvmrname='kp20'//specname(1:6)//'.vmr'
         endif
      elseif(specname(1:2).eq.'ss' .or. specname(1:2).eq.'sr') then
         if(lspn.lt.6) stop 'lspn<6'
         newmodname=specname(1:6)//'.mod'    ! ACE
         newvmrname=specname(1:6)//'.vmr'    ! ACE
      elseif(specname(1:2).eq.'s1' .or. specname(1:2).eq.'r1') then
         if(lspn.lt.6) stop 'lspn<6'
         newmodname=specname(1:6)//'.mod'    ! ACE
         newvmrname=specname(1:6)//'.vmr'    ! ACE
      elseif(specname(1:24).eq.'FTS-TRANSMITTANCE-NA-NA-') then
         newmodname=rlgfile(1:lr-4)//'.mod'    ! ACE (new names)
         newvmrname=rlgfile(1:lr-4)//'.vmr'    ! ACE (new names)
      elseif(rlgfile(lr-1:lr).eq.'ws') then
         if(lspn.lt.6) stop 'lspn<6'
         newmodname=specname(2:6)//'.zpt' ! Wollongong spectra
         newvmrname=specname(2:6)//'.ref' ! Wollongong spectra
      elseif(specname(1:2).eq.'iz') then
         if(lspn.lt.6) stop 'lspn<6'
         newmodname='pt_'//specname(1:10)//'.prf' ! Izana
         newvmrname='vmr_'//specname(1:10)//'.ref' ! Izana
      elseif(specname(1:2).eq.'o2') then
         if(lspn.lt.6) stop 'lspn<6'
         newmodname=specname(1:6)//'av.mod' ! TMF FTUVS O2 A-band spectra
         newvmrname=specname(1:6)//'av.vmr' ! TMF FTUVS O2 A-band spectra
      elseif(specname(1:1).eq.'9' .or. specname(1:1).eq.'0') then
         if(lspn.lt.6) stop 'lspn<6'
         newmodname=specname(1:6)//'.mod'
         newvmrname=specname(1:6)//'.vmr'
      elseif(specname(1:3).eq.'j60' ) then
         newmodname=specname(1:lspn)//'.mod'
         newvmrname=specname(1:lspn)//'.vmr'
      else
         if(lspn.lt.10) stop 'lspn<10'
         newmodname=specname(1:10)//'.mod' !  CIT
         newvmrname=specname(1:10)//'.vmr' !  CIT
      endif
      lm=lnbc(newmodname)
      lv=lnbc(newvmrname)

      if(newmodname(:lm).eq.modname(:lm)) then  ! Current model okay
         newmod=.false.
      else                       ! Look for new model
c         write(6,*) ' Searching for ZPT model ',newmodname
         inquire(file=gggdir(:lrt)//'models'//dl//ext//dl//newmodname,
     $   exist=newmod)
         lmn=lnbc(modname)
         if(newmod) then
            modname=newmodname
            lmn=lnbc(modname)
            write(lunw_rpt,*)'Reading ',modname(:lmn),' for ',specname

c           read(modname,'(2x,i4,i2,i2)') iyyy,im,id
c           date_mod=iyyy+float(im-1)/12.0+id/365.25
c           iday=nint(30.5*float(im-1)+id)
c           write(*,*)'date_mod=',date_mod
            date_mod=iyr+idoy/365.25 !from the runlog
            iday=idoy !from the runlog
c           write(*,*)'Setting date_mod to runlog value: ',date_mod

            call read_model_fc(lunr_mod,
     &      gggdir(:lrt)//'models'//dl//ext//dl//modname,
c     $      z,t,p,d,h2ovmr,mmw,roc,nlev,ztrop_ncep,ztrop_gct)
     $      z,mmw,nlev,t,p,d,h2ovmr,roc,ztrop_ncep,ztrop_gct,tlat)
         elseif(lmn.eq.0) then        !first spectrum in runlog
            write(lunw_rpt,*)'Unable to find ZPT model '//newmodname
            write(*,*)'Unable to find ZPT model '//newmodname
            menuinput=gggdir(:lrt)//'models'//dl//ext//dl//'models.men'
            call readmenu(lunr_men,menuinput,modname)
            write(lunw_rpt,*)'Reading ',modname(:lmn),' for ',specname
            call read_model_fc(lunr_mod,
     &      gggdir(:lrt)//'models'//dl//ext//dl//modname,
c     $      z,t,p,d,h2ovmr,mmw,roc,nlev,ztrop_ncep,ztrop_gct)
     $      z,mmw,nlev,t,p,d,h2ovmr,roc,ztrop_ncep,ztrop_gct,tlat)
            newmod=.true.
            date_mod=iyr+idoy/365.25 ! Set to runlog value 
c           write(*,*)'Setting date_mod to runlog value: ',date_mod
         else
            write(lunw_rpt,*) ' Cannot find ZPT model '//newmodname//
     &      '  Will re-use: '//modname
         endif
      endif  !  (newmodname.ne.modname)
c
      if(newvmrname(:lv).eq.vmrname(:lv)) then ! Existing vmr okay
         newvmr=.false.
      else   ! search for new vmr profiles
c         write(6,*) ' Searching for vmr profiles ',newvmrname
         inquire(file=gggdir(:lrt)//'vmrs'//dl//ext//dl//newvmrname,
     $   exist=newvmr)
         lvn=lnbc(vmrname)
         if(newvmr) then
            vmrname=newvmrname
            lvn=lnbc(vmrname)
            write(lunw_rpt,*)'Reading ',vmrname(:lvn),' for ',specname
            vmr_found=.true.
         else if(lvn.eq.0) then ! First spectrum
            write(lunw_rpt,*)
     &            ' Unable to find vmr profiles '//newvmrname
            write(*,*) ' Unable to find vmr profiles '//newvmrname
            call readmenu(lunr_men,
     &      gggdir(:lrt)//'vmrs'//dl//ext//dl//'vmrs.men',vmrname)
          write(lunw_rpt,*)'Reading ',vmrname(:lvn),' for ',specname
            newvmr=.true.
            vmr_found=.false.
         else
            write(lunw_rpt,*) ' Cannot find vmr profile '//
     &      newvmrname(:lnbc(newvmrname))//
     &      '  Will re-use: '//vmrname(:lvn)
            vmr_found=.false.
         endif
c         call read_vmrs_fc(lun_vmr,
c     &   root(:lrt)//'vmrs'//dl//ext//dl//vmrname,z,nlev,vmrlabel,
c     &  ztrop_gct,apvmr,mgas,root(:lrt)//'models'//dl//ext//dl//modname,
c     &  tlat,idoy,vmr_found)
        call read_refvmrs(lunr_vmr,
     &   gggdir(:lrt)//'vmrs'//dl//ext//dl//vmrname,nlev,z,mgas,
     &   gggdir(:lrt)//'models'//dl//ext//dl//modname,
     &   vmrlabel,refvmr,ngas,reflat_vmr,date_vmr,ztrop_vmr)

c If there's only 1 level (i.e. lab), or the vmr file for that particular
c measurement already exists, don't try to modify the .vmr file
        if ((nlev.gt.1).and.(vmr_found.eqv..false.)) then
        call resample_vmrs_at_effective_altitudes(nlev,z,mgas,ngas,
     &   refvmr,ztrop_gct,ztrop_vmr,oblat,reflat_vmr,apvmr)
         call apply_vmr_latitude_gradients(nlev,z,mgas,apvmr,
     &   ztrop_gct,reflat_vmr,oblat,apvmr)
         call apply_secular_trends(nlev,z,mgas,apvmr,
     &   ztrop_gct,reflat_vmr,oblat,date_mod,date_vmr,apvmr)
         fryr=date_mod-int(date_mod)
         do ilev=1,nlev
            z8=dble(z(ilev))
            do jgas=1,6
            apvmr(jgas,ilev)=apvmr(jgas,ilev)*
     &      compute_seasonal_cycle(jgas,z8,ztrop_gct,oblat,fryr)
            end do
         end do
        else
           apvmr=refvmr
        endif

      endif   !  if(newvmrname.ne.vmrname) then
c
c  Output model information (SUNRUN.MAV)
      if(newvmr.or.newmod) then
         t_cell(1)=tins+273.16
         write(lunw_mav,'(a)')'Next Spectrum:'//specname
         write(lunw_mav,'(3i4)') 3, nspeci+4, nlev+ncell
         write(lunw_mav,'(a,2f9.3)') 'Tropopause Altitudes (NCEP/GCT):',
     &    ztrop_ncep,ztrop_gct
         zpbl=0.0
         if (vmr_found) then
c             Do nothing when the vmr file is found.
c             write(*,*)'vmr_found, ',vmr_found,vmrname(:lvn)
c             write(71,*)apvmr
         else
c  For ground-based observations, over-write apriori H2O vmr with NCEP model.
c  For HDO, the factor 0.15*(8.0+log10(h2ovmr(ilev))) adjusts for the effects
c  of isotopic fractionation.
c  For H2O = 0.03,  HDO=0.975*H2O (sea-level in tropics)
c  For H2O = 0.01,  HDO=0.900*H2O (sea-level at mid-latitude)
c  For H2O = 3E-06, HDO=0.375*H2O (tropopause)
         if(h2ovmr(1).gt.0.0) then ! An H2O profile was found in the .mod file
            if ( pout.gt.600.0 ) then ! ground-based observations
c               write(*,*) 'Replacing H2O & HDO vmrs with NCEP profiles'
               do ilev=1,nlev
                 apvmr(1,ilev)=h2ovmr(ilev)  ! H2O
                 apvmr(49,ilev)=h2ovmr(ilev)*
     &             0.16*(8.0+log10(h2ovmr(ilev))) ! HDO
               end do
            endif
         endif

c              call setvmr(apvmr,mgas,z,h2ovmr,co2vmr,nlev,iyr,idoy,
c     &        zpdtim,oblat,oblon,obalt,tout,pout,hout,ztrop_gct,zpbl)
         endif
         call write_mav(z,t,p,d,apvmr,nlev,lunw_mav,ncell,
     &   t_cell,p_cell,vmr_cell,gas_in_cell,vmrlabel,isofile,mgas)
      endif
c------------------------------------------------------------------
      if(pout.eq.0.0) then
         zpres=obalt
      else
         zpres=height(sngl(pout)/1013.25,p,t,z,nlev)
      endif
c
      if(abs(obalt-zpres).gt.0.1) write(lunw_rpt,*)
     $'Warning: Geometrical & Pressure altitudes differ: ',obalt,zpres
c
      do ilev=1,nlev
        r(ilev)=z(ilev)+roc
      end do
c  Determine ALT from POBS if > 400.0    Otherwise use geometric altitude
      if(pout.lt.400.0d0 .and. iyr.lt.1979) then
        zobs=obalt
      else
        zobs=zpres   ! prefer to determine GRADIUS from POBS
      endif
c
c  Determine ray slant paths
      fovr=90.*sngl(fovo)/spi     ! convert radians diameter to deg radius
c      write(*,*)'Calling TLPATH....',solzen,zobs,specname(:20)
      call tlpath(nlev,z,t,p,solzen,fovr,roc,zobs,
     $ sngl(wavtkr),frqcen,zmin,bend,splos,rc)
      if(rc.ne.0) write(6,*)'Error in TLPATH:',rc,specname,solzen
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
c         if(nlev.gt.1) stop 'vmrs dont extend high enough'
         splos(nlev)=splos(nlev)+aipl*(pout/1013.25)/p(nlev)
         if(nlev.eq.1) splos(nlev)=aipl
       else
         delta_p=p(ilev-1)-p(ilev)
         splos(ilev-1)=splos(ilev-1)+aipl*(pout/1013.25-p(ilev))/delta_p
         splos(ilev)=splos(ilev)+aipl*(p(ilev-1)-pout/1013.25)/delta_p
       endif
c
      write(lunw_ray,'(a57,260(1x,f10.5))')specname,obalt,pout,
     & solzen,bend,fovo,zmin,(cell_length(jcell),jcell=1,ncell),
     & (wlimit(dble(multipath)*splos(j),'f10.5'),j=1,nlev)
      write(lunw_rpt,'(a57,16f11.5)')specname,obalt,pout,solzen,
     & bend,fovo,zmin,(cell_length(jcell),jcell=1,ncell) 
      go to 1
 99   close(lunr_rlg)
      close(lunw_ray)
      close(lunw_mav)
      close(lunw_rpt)
      if(nspe.le.0)write(6,*)'none of these spectra could be accessed'
      if(nspe.le.0)write(6,*)'or error reading runlog (wrong format?)'
      write(*,*) nspe,' spectra total'

c  Code to generate post_processing.sh batch file and associated inputs.

      call write_postprocessfile(ext,rlgfile)

      stop
      end
