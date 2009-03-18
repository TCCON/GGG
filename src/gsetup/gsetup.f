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
      integer*4 apo,bytepw,fbc,fnbc,ifirst,ilast,ilev,
     & interp,iset,iyr,j,k,w1,w2,le,lg,lc,lnbc,lrt,lunr,lun_ray,lun_ggg,
     & lun_rlg,lun_mul,lun_pp,lun_mav,lun_rpt,lun_mod,lun_vmr,lunz,
     & mlev,nlev,nlhead_rlg,nlhead_ray,nlhead_ggg,
     & lr,la,mgas,i,ncol,
c    & platform,
     & multipath, nspe,possp,rc,istat,lm,lv,gas_in_cell,
     & lmn,lvn,lun_men
      parameter (lun_men=12)     ! for reading the various .men files
      parameter (lun_ray=13)     ! for writing the .ray file
      parameter (lunr=14)        ! for reading (general purpose) 
      parameter (lun_ggg=15)     ! for writing the .ggg files
      parameter (lun_mul=16)     ! for writing the multiggg.sh file
      parameter (lun_mav=17)     ! for writing the .mav file
      parameter (lun_rpt=18)     ! for writing gsetup.rpt
      parameter (lun_rlg=21)     ! for reading the runlog
      parameter (lun_mod=22)     ! for reading the .mod file (readmodFC)
      parameter (lun_vmr=23)     ! for reading the .vmr file (readvmrFC)
      parameter (lunz=24)        ! for reading (general purpose) 
      parameter (lun_pp=25)      ! for post_processing.sh
      parameter (mlev=250)       ! maximum number of atmospheric levels
      parameter (mgas=80)        ! maximum number of gases
      parameter (nlhead_ray=3)    ! Number of header lines in .ray file.
      parameter (nlhead_ggg=18)  ! Number of header lines in .ggg file.
      real*4 erroff
      parameter (erroff=0.004,apo=2,interp=1)
      logical*4 newmod,newvmr

      real*4 dum,freq,frqcen,height,
     & zmin,zpres,p_cell,t_cell,
     & zobs,             ! the geocentric radius of the observer fed to TLPATH
     & bend,             ! total ray bending due to refraction
     & d(mlev),          ! number densities at levels (molec.cm-3)
     & h2ovmr(mlev),     ! h2o profile from model file
     & co2vmr(mlev),     ! co2 profile from simulate_co2_vmr
     & fovr,             ! External angular radius of FOV in degrees
     & mmw(mlev),        ! mean molecular weights
     & p(mlev),          ! pressures of levels (atm.)
     & r(mlev),          ! radii of levels (km) = z + roc
     & roc,              ! radius of curvature (km)
     & cell_length,      ! Length of cell permanently in solar beam
     & solzen,           ! solar zenith angle used (may be refracted)
     & splos(mlev),      ! array of line-of-sight slant paths
     & t(mlev),          ! temperatures of levels (K)
     & vmr(mgas,mlev),   ! buffer for vmr's
     & z(mlev)           ! altitudes of levels (km)

      real*8  ppbl_atm, ptrop, ptrop_atm, delta_p, wlimit,
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
     & fvsi,             ! Fractional Variation in Solar Intensity
     & wspd,             ! Wind Speed (m/s)
     & wdir,             ! Wind Direction (deg)
     & aipl,             ! Airmass-Independent Path length (km)
     & lasf,             ! Laser Frequency (e.g. 15798 cm-1)
     & wavtkr            ! Frequency of the sun-tracker (actively-tracked)

        
      character
     & apf*2,              ! apodization function (e.g. 'BX','TR' etc)
     & dl*1,               ! delimiter (='/' Unix, ='\' DOS)
     & og*1,               ! observation geometry
     & ext*3,              ! observation geometry
     & filnam*48,          ! general purpose character buffer
     & filnamwas*48,       ! general purpose character buffer
     & header*92,          ! general purpose character buffer
     & levels*24,          ! name of file containing levels
     & listof*24,          ! the selected list of windows
     & modname*80,         ! path to new model
     & newmodname*80,      ! path to new model
     & vmrname*80,         ! path to new vmr
     & newvmrname*80,      ! path to new vmr
     & menuinput*80,       ! path to xxx.men file
     & prvwin*14,          ! name of previous window (i.e. GAS_1234)
     & root*64,            ! root directory
     & specname*35,        ! name of spectrum
     & runlog*40,          ! name of occultation file
     & slk*1,              ! symbol which links gas name and frequency
     & col1,               ! 
     & ccell*1,            ! Character in specname denoting internal cell 
c     & user*8,             ! investigator
     & version*46,         ! program version number
     & isofile*80,         ! path to "isotopolog.dat"
     & vmrlabel*1000,      ! column labels from vmr file
     & window*120          ! name of window

      version= ' GSETUP Version 2.7.5       6 Mar 2009    GCT '
      modname='                                                '
      vmrname='                                                '
      filnamwas='qwertyuioqwertyuioqwertyuioqwertyuio'
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cc     Platform specification:                                                    DG000909
c      call getenv('LOGNAME',user)
c      if(user.eq.'        ')then ! It's not a Unix/Linux machine
c         platform=2               !2=PC-Win32
c         dl=char(92)              !Back-slash  ('\')
c         root='g:'//dl            !DG Jan03
c         user='PC-Win'
c      else                       ! It's a Unix/Linux machine
c         platform=0            
c         dl='/'
c         call getenv('GGGPATH',root)
c         root=root(:lnbc(root))//dl
cc         root='/home/toon/ggg/'
c      endif
c      lrt=lnbc(root)       !Length of root
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      dl='/'
      call getenv('GGGPATH',root)
      root=root(:lnbc(root))//dl
      lrt=lnbc(root)       !Length of root
c---------------------------------------------------------------------
      slk='_'
      ext(1:3)='   '
      open(lun_rpt,file='gsetup.rpt', status='unknown')

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

      if(ext(1:1).eq.'b' .or. ext(1:1).eq.'b') then
          og='l'
      else
          og='v'
      endif
c------------------------------------------------------------------
c  Choose which runlog to analyze
      menuinput=root(:lrt)//'runlogs'//dl//ext//dl//'runlogs.men'
      call readmenu(lun_men,menuinput,runlog)
      lr=lnbc(runlog)
c------------------------------------------------------------------
c  choose which levels to use
      menuinput=root(:lrt)//'levels'//dl//'levels.men'
      call readmenu(lun_men,menuinput,levels)
c  Read chosen level altitudes
      open(lunr,file=root(:lrt)//'levels'//dl//levels(:lnbc(levels)),
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
      menuinput=root(:lrt)//'windows'//dl//ext//dl//'windows.men'
      call readmenu(lun_men,menuinput,listof)

c  read the individual microwindows and create the xxxxxxxx.ggg input files
      open(lun_mul,file='multiggg.sh',status='unknown')     
      open(lunr,file=root(:lrt)//'windows'//dl//ext
     $ //dl//listof(:lnbc(listof)),status='old')
      read(lunr,'(a91)') header
      prvwin='            '
      do k=1,999999
 666    read(lunr,'(a)',end=103) window
        if(window(1:1).eq.':' .or. window(1:1).eq.';') go to 666
        w1=fnbc(window)
        if(w1.eq.0) go to 666 ! blank line
        w2=index(window,'.')
        read(window(w1:w2),*) freq
        write(6,'(a)') ' '//window(:lnbc(window))

        lc=index(window,':')  ! lc=position of colon
        lg=lc+fnbc(window(lc+1:))
        le=lg+fbc(window(lg+1:))-1
        filnam=window(lg:le)//slk//window(w1:w2)//runlog(:lr-3)//'ggg'
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
        open(lun_ggg,file=filnam,status='unknown')
        write(lun_ggg,'(i3)') nlhead_ggg
        write(lun_ggg,'(a)') version
c        write(lun_ggg,'(a)') version
        write(lun_ggg,'(a)') root(:lrt)//'config'//dl//'data_part.lst'
        write(lun_ggg,'(a)') root(:lrt)//'apriori'//dl//'gfit_ap.'//ext
        write(lun_ggg,'(a)') root(:lrt)//'runlogs'//dl//ext//dl//
     $  runlog(:lr)
        write(lun_ggg,'(a)') root(:lrt)//'levels'//dl//levels
        write(lun_ggg,'(a)') root(:lrt)//'models'//dl//ext//dl
        write(lun_ggg,'(a)') root(:lrt)//'vmrs'//dl//ext//dl
        write(lun_ggg,'(a)') runlog(:lr-3)//'mav'
        write(lun_ggg,'(a)') runlog(:lr-3)//'ray'
        write(lun_ggg,'(a)') root(:lrt)//'isotopologs'//dl//
     &  'isotopologs.dat'
        write(lun_ggg,'(a)') root(:lrt)//'windows'//dl//ext//dl//listof

c  write names of linelists
c        write(lun_ggg,'(a)')      ! llsize
        write(lun_ggg,'(a)')
     &  root(:lrt)//'linelist'//dl//'atm.101 '//
     &  root(:lrt)//'linelist'//dl//'gct.101 '//
     &  root(:lrt)//'linelist'//dl//'fcia.101 '//
     &  root(:lrt)//'linelist'//dl//'scia.101'
        write(lun_ggg,'(a)') root(:lrt)//'linelist'//dl//'solar_dc.101'

c  write location of .ak files (averaging kernels)
        write(lun_ggg,'(a)') root(:lrt)//'ak'//dl//'k'

c  write location of spectral fits
        write(lun_ggg,'(a)') root(:lrt)//'spt'//dl//'z'

        write(lun_ggg,'(a)') filnam(:lnbc(filnam)-3)//'col'
c        write(lun_ggg,'(a)') filnam(lnbc(filnam)-11:lnbc(filnam)-3)//'col'
        write(lun_ggg,'(a)') window
        close(lun_ggg)
c
        write(lun_mul,'(a)') root(:lrt)//'bin'//dl//'gfit<'//
     $  filnam(:lnbc(filnam))//'>/dev/null'
      end do   !  k=1,999999
 103  close(lunr)
      close(lun_mul)
c------------------------------------------------------------------
c  Compute the slant paths and write them to disk.
      write(lun_rpt,*)
      write(lun_rpt,'(a)')
     $'  Spectrum        Zobs   Pobs    ASZA    BEND    FOVO    TANG'
      write(lun_rpt,'(a)')
     $'                   km    mbar    deg.    deg.    mrad     km'
      nspe=0
      open(lun_mav,file=runlog(:lr-3)//'mav',status='unknown')
c  Open the file which will contain the slant paths
      write(lun_mav,'(a)')version
      write(*,*)' Opening: ',runlog(:lr-3)//'ray'
      open(lun_ray,file=runlog(:lr-3)//'ray',status='unknown')
      write(lun_ray,'(i2,i4)') nlhead_ray,nlev+8
      write(lun_ray,'(a)') version
c      write(lun_ray,'(a38,<mlev>(a2,i3.3))')
      write(lun_ray,'(a,250(a2,i3.3))')
     $' SpectrumName      Zobs  Pobs  ASZA  Bend  FOV  Zmin  Cell',
     & (' L',j,j=0,nlev-1)
c  Get the header/runlog information pertaining to RUNLAB
      write(6,*) ' Runlog=',root(:lrt)//'runlogs'//dl//ext//dl//
     & runlog(:lr)
      open(lun_rlg,file=root(:lrt)//'runlogs'//dl//ext//dl//
     & runlog,status='old')
      read(lun_rlg,*) nlhead_rlg,ncol           ! read runlog header line
      do i=2,nlhead_rlg
         read(lun_rlg,*)
      end do
1     call read_runlog(lun_rlg,col1,specname,iyr,iset,zpdtim,oblat,
     & oblon,obalt,asza,zenoff,azim,osds,
     & opd,fovi,fovo,amal,ifirst,ilast,
     & graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     & tout,pout,hout,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)
       la=lnbc(specname)
c       write(*,*)istat,la,' '//specname
       if(istat.ne.0) go to 99   !  EOF
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
      if(mod(nspe,1000).eq.0) write(*,*) nspe,' spectra'

c  Figure out whether the spectrum was through an internal cell.
      p_cell=0.0
      cell_length=0.0
      if(la.ge.20) then
         ccell=specname(12:12)
c         if( ccell.eq.'a' .or. ccell.eq.'b' .or. ccell.eq.'c') then
         if( ccell.ne.'0') then
            cell_length=0.0001 ! km (10cm)
            p_cell=5.097       ! mbar ! OCO Brukers have internal cells
            gas_in_cell=15     ! HCl
         endif
      endif
c
      if(specname(2:3).eq.'hg' .or. specname(2:3).eq.'in') then
         if(la.lt.8) stop 'la<8'
c         newmodname=runlog(1:5)//specname(4:8)//'.mod'  ! MkIV
c         newvmrname=runlog(1:5)//specname(4:8)//'.vmr'  ! MkIV
         newmodname=ext//specname(4:8)//'.mod'  ! MkIV
         newvmrname=ext//specname(4:8)//'.vmr'  ! MkIV
      elseif(specname(7:9).eq.'R0.') then  ! Kitt Peak
         if(specname(1:1).eq.'7' .or. specname(1:1).eq.'8' .or.
     &   specname(1:1).eq.'9' ) then
         if(la.lt.8) stop 'la<8'
         newmodname=specname(1:8)//'.mod'
         newvmrname=specname(1:8)//'.vmr'
         endif
      elseif(specname(1:2).eq.'ss' .or. specname(1:2).eq.'sr') then
         if(la.lt.6) stop 'la<6'
         newmodname=specname(1:6)//'.mod'    ! ACE
         newvmrname=specname(1:6)//'.vmr'    ! ACE
      elseif(specname(1:2).eq.'s1' .or. specname(1:2).eq.'r1') then
         if(la.lt.6) stop 'la<6'
         newmodname=specname(1:6)//'.mod'    ! ACE
         newvmrname=specname(1:6)//'.vmr'    ! ACE
      elseif(runlog(lr-1:lr).eq.'ws') then
         if(la.lt.6) stop 'la<6'
         newmodname=specname(2:6)//'.zpt' ! Wollongong spectra
         newvmrname=specname(2:6)//'.ref' ! Wollongong spectra
      elseif(specname(1:2).eq.'iz') then
         if(la.lt.6) stop 'la<6'
         newmodname='pt_'//specname(1:10)//'.prf' ! Izana
         newvmrname='vmr_'//specname(1:10)//'.ref' ! Izana
      elseif(specname(1:2).eq.'o2') then
         if(la.lt.6) stop 'la<6'
         newmodname=specname(1:6)//'av.mod' ! TMF FTUVS O2 A-band spectra
         newvmrname=specname(1:6)//'av.vmr' ! TMF FTUVS O2 A-band spectra
      elseif(specname(1:1).eq.'9' .or. specname(1:1).eq.'0') then
         if(la.lt.6) stop 'la<6'
         newmodname=specname(1:6)//'.mod'
         newvmrname=specname(1:6)//'.vmr'
      else
         if(la.lt.10) stop 'la<10'
         newmodname=specname(1:10)//'.mod' !  CIT
         newvmrname=specname(1:10)//'.vmr' !  CIT
      endif
      lm=lnbc(newmodname)
      lv=lnbc(newvmrname)

      if(newmodname(:lm).eq.modname(:lm)) then  ! Current model okay
         newmod=.false.
      else                       ! Look for new model
c         write(6,*) ' Searching for ZPT model ',newmodname
         inquire(file=root(:lrt)//'models'//dl//ext//dl//newmodname,
     $   exist=newmod)
         lmn=lnbc(modname)
         if(newmod) then
            modname=newmodname
            lmn=lnbc(modname)
            write(lun_rpt,*)'Reading ',modname(:lmn),' for ',specname
            call readmodFC(lun_mod,
     &      root(:lrt)//'models'//dl//ext//dl//modname,
     $      z,t,p,d,h2ovmr,mmw,roc,nlev,ptrop)
         elseif(lmn.eq.0) then        !first spectrum in runlog
            write(lun_rpt,*)'Unable to find ZPT model ',newmodname
            write(*,*)'Unable to find ZPT model ',newmodname
            menuinput=root(:lrt)//'models'//dl//ext//dl//'models.men'
            call readmenu(lun_men,menuinput,modname)
            write(lun_rpt,*)'Reading ',modname(:lmn),' for ',specname
            call readmodFC(lun_mod,
     &      root(:lrt)//'models'//dl//ext//dl//modname,
     $      z,t,p,d,h2ovmr,mmw,roc,nlev,ptrop)
            newmod=.true.
         else
            write(lun_rpt,*) ' Cannot find ZPT model '//newmodname//
     &      '  Will re-use: '//modname
         endif
      endif  !  (newmodname.ne.modname)
c
      if(newvmrname(:lv).eq.vmrname(:lv)) then ! Existing vmr okay
         newvmr=.false.
      else   ! search for new vmr profiles
c         write(6,*) ' Searching for vmr profiles ',newvmrname
         inquire(file=root(:lrt)//'vmrs'//dl//ext//dl//newvmrname,
     $   exist=newvmr)
         lvn=lnbc(vmrname)
         if(newvmr) then
            vmrname=newvmrname
            lvn=lnbc(vmrname)
            write(lun_rpt,*)'Reading ',vmrname(:lvn),' for ',specname
            call readvmrFC(lun_vmr,
     &      root(:lrt)//'vmrs'//dl//ext//dl//vmrname,z,nlev,vmrlabel,
     &      vmr,mgas,root(:lrt)//'models'//dl//ext//dl//modname)
         else if(lvn.eq.0) then ! First spectrum
            write(lun_rpt,*) ' Unable to find vmr profiles ',newvmrname
            write(*,*) ' Unable to find vmr profiles ',newvmrname
            call readmenu(lun_men,
     &      root(:lrt)//'vmrs'//dl//ext//dl//'vmrs.men',vmrname)
          write(lun_rpt,*)'Reading ',vmrname(:lvn),' for ',specname
            call readvmrFC(lun_vmr,
     &      root(:lrt)//'vmrs'//dl//ext//dl//vmrname,z,nlev,vmrlabel,
     &      vmr,mgas,root(:lrt)//'models'//dl//ext//dl//modname)
            newvmr=.true.
         else
            write(lun_rpt,*) ' Cannot find vmr profile '//
     &      newvmrname(:lnbc(newvmrname))//
     &      '  Will re-use: '//vmrname(:lvn)
         endif
      endif   !  if(newvmrname.ne.vmrname) then
c
c  Output model information (SUNRUN.MAV)
      isofile=root(:lrt)//'isotopologs/isotopologs.dat'
c      write(*,*) specname,newmod,newvmr
      if(newvmr.or.newmod) then
         t_cell=tins+273.16
         write(lun_mav,'(a)')'Next Spectrum:'//specname
         ppbl_atm=0.0
         ptrop_atm=ptrop/1013.25
c         write(*,*)modname,h2ovmr(1),h2ovmr(2),h2ovmr(3),h2ovmr(4)
         call setvmr(vmr,mgas,z,t,p,h2ovmr,co2vmr,nlev,iyr,iset,
     &   zpdtim,oblat,oblon,obalt,tout,pout,hout,ptrop_atm,ppbl_atm)
         call write_mav(z,t,p,d,vmr,nlev,lun_mav,
     &   t_cell,p_cell,gas_in_cell,vmrlabel,isofile,mgas)
      endif
c------------------------------------------------------------------
      if(pout.eq.0.0) then
         zpres=obalt
      else
         zpres=height(sngl(pout)/1013.25,p,t,z,nlev)
      endif
c
      if(abs(obalt-zpres).gt.0.1) write(lun_rpt,*)
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
      fovr=90.*sngl(fovo)/3.14159265 ! convert radians diameter to deg radius
c      write(*,*)'Calling TLPATH....',solzen,zobs
      call tlpath(nlev,z,t,p,solzen,fovr,roc,zobs,
     $sngl(wavtkr),frqcen,zmin,bend,splos,rc)
      if(rc.ne.0) write(6,*)'Error in TLPATH:',rc
c
c  Distribute AIPL between the two atmospheric levels that bracket POUT
       do ilev=2,nlev
          if(pout/1013.25.gt.p(ilev)) exit
       end do
       if(ilev.gt.nlev) then  ! Failed to find low enough pressure level
         if(aipl.gt.0) then
            write(*,*) 'Warning: aipl>0; pout<p(i)'
            write(*,*) aipl,pout,p(nlev),nlev
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
      write(lun_ray,'(a35,7f10.4,250(1x,f9.4))')specname,obalt,pout,
     & solzen,bend,fovo,zmin,cell_length,
     & (wlimit(dble(multipath)*splos(j),'f9.4'),j=1,nlev)
      write(lun_rpt,'(a35,7f10.4)')specname,obalt,pout,solzen,
     & bend,fovo,zmin,cell_length 
      go to 1
 99   close(lun_rlg)
      close(lun_ray)
      close(lun_mav)
      close(lun_rpt)
      if(nspe.le.0)write(6,*)'none of these spectra could be accessed'
      if(nspe.le.0)write(6,*)'or error reading runlog (wrong format?)'
      write(*,*) nspe,' spectra'

c  Code to generate post_processing.sh batch file and associated inputs.

      open(lun_pp,file='post_processing.sh',status='unknown')     

      write(lun_pp,'(a)')'~/ggg/bin/collate_results<.collate.input'
      open(lunz,file='.collate.input', status='unknown')
      write(lunz,'(a)') og
      close(lunz)

      write(lun_pp,'(a)')'~/ggg/bin/average_results<'//
     &'.average_results.input'
      open(lunz,file='.average_results.input',status='unknown')
      write(lunz,'(a)') runlog(:lr-3)//og//'sw'
      close(lunz)

      write(lun_pp,'(a)') '~/ggg/bin/apply_airmass_correction<'//
     &'.apply_airmass_correction.input'
      open(lunz,file='.apply_airmass_correction.input',status='unknown')
      write(lunz,'(a)') runlog(:lr-3)//og//'av'
      close(lunz)

      write(lun_pp,'(a)') '~/ggg/bin/apply_insitu_correction<'//
     &'.apply_insitu_correction.input'
      open(lunz,file='.apply_insitu_correction.input',status='unknown')
      write(lunz,'(a)') runlog(:lr-3)//og//'av.ada'
      close(lunz)

      write(lun_pp,'(a)') '~/ggg/bin/write_official_output_file<'//
     &'.write_official_output_file.input'
      open(lunz,file='.write_official_output_file.input',
     & status='unknown')
      write(lunz,'(a)') runlog(:lr-3)//og//'av.ada.aia'
      close(lunz)

      close(lun_pp)
      stop
      end
