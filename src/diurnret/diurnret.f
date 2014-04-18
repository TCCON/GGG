c  Program diurnret
c  Retrieves VMR profiles from slant column abundances.
c
c  Input Files:
c      runlog.lav       measured slant columns abundances,
c      runlog.ray       matrix ofcalculated matrix of slant path distances.
c      balYYDDD.mav     T, P, D, VMR profiles for each date
c      dcfYYDDD.gas     Diurnal correction matrix (optional)
c
c  Output files:
c      runlog.lav.ret   retrieved vmr profiles for each date
c
c  Version History:
c  ' V1.0.0    Oct-92  GCT' Program used for analysis of Sep 92 balloon flight.
c  ' V3.0.0 07-Dec-95  GCT'
c  ' V3.0.1 23-Mar-96   BS' increased mspeci from 62 to 65
c  ' V3.1.0 08-Apr-96  GCT' changed read format of .lav file
c  ' V4.0.0 13-Apr-96  GCT' included diurnal corrections
c  ' V4.1.1 26-Apr-96   BS' increased msza to 13
c  ' V4.2.0 19-Mar-97  GCT' Made SMOO information-dependent
c  ' V4.3.0 24-Mar-97  GCT' Interpolated UT,SOLZEN,TLAT,TLON to
c                           retrieval altitudes
c  ' V4.4.0 31-Mar-97  GCT' Supports double column labels
c  ' V4.4.1 22-Jun-97  GCT' Removed () from column labels
c  ' V4.4.2 23-Jun-97  GCT' Increased MSZA to 17
c  ' V4.4.3 24-Jul-97  GCT' New version of SUBSTR.F
c  ' V4.4.4 15-Sep-97  GCT' Improved interpolation of TLAT,TLON etc
c  ' V4.4.5 10-Oct-97  GCT' New format for DCF, VMR, & MAV files
c  ' V4.4.6  1-Feb-98  GCT' Corrected format statement to
c  74     format(a14,5f10.4,150f8.3)   from     format(a14,5f9.4,150f7.3)
c  ' V4.5.0  2-Feb-98  GCT' Uses ZENAZ instead of SUBSOLAR
c  ' V4.5.1  3-Feb-98  GCT' Fixed error in NCOL in .ret file
c  ' V4.5.2  4-Feb-98  GCT' Uses .in instead of .hg SP file.
c  ' V4.5.3  5-Feb-98  GCT' Increased PABEL & STRING to C*9600
c  ' V4.5.4  6-Feb-98  GCT' Output file is OCCNAME.lxx.ret (not OCCNAME.ret)
c  ' V4.6.0  7-Feb-98  GCT' OCCNAME.ATS output file no longer used.
c  ' V4.6.1 12-Feb-98  GCT' reduced smoothing from /12 to /16
c  ' V4.6.2 20-Feb-98  GCT' now recognizes "no_1903" and "no-1903"
c  ' V4.6.3 13-May-98   BS' inserts date+time and retrieval range
c  ' V4.7.0 20-Sep-98  GCT' New format (YEAR first) .lav and .ret files.
c  ' V4.8.0  6-Oct-98  GCT' Uses effective altitude (ZEFF) in determinations...
c  ' V4.9.0 14-Oct-98  GCT' Copies previous headers (from GFIT, GGGAVG).
c  ' V4.9.1 21-Jan-99  GCT' Increased dimension and length of COMMENT.
c  ' V4.9.2 15-May-99  GCT' Fixed bug in ZENAZ which produced wrong TPLAT & TPLONG
c  ' V4.9.3 18-May-99  GCT' New scheme for int/extrapolating to retrieval altitude
c  ' V4.9.4  7-Sep-99  GCT' Increased mspeci from 65 to 122 for isotopic retieval
c  ' V4.9.5  8-Sep-99  GCT' Uses column labels from .mav file (instead of MOLNUM) to identify gases.
c  ' V4.9.7  5-Jan-00  GCT' Increased mspeci from 100 to 122 for isotopic retieval
c  ' V4.9.8 28-May-00  GCT' Defined LTOP variable to simplify code.
c  ' V4.9.9 19-Jul-00  GCT' Added IF's to open right input files for ATMOS
c  ' V4.9.9  9-Jan-01  GCT' Increased MCOMM to 6
c    V4.9.10  29-Oct-2001   GCT' Added "d0" to DMOD statements. 9000 format
c    V4.10.1  29-Jul-2002   Calculate Minimum Detectable Abundance vmr (Vmda)
c    V4.10.2  30-Jul-2002   GCT' Added "diurnret.out" output file.
c    V4.10.3  23-Jul-2003   GCT' Increased parameter msza from 22 to 130
c    V4.10.3  23-Jul-2003   Increased dimension of dcflabel from 256 to 1024
c    V4.10.4   1-Nov-2004   Increased parameter mspeci from 126 to 140
c    V4.10.5   4-Mar-2005   Deleted  calls to clistindex (already commented out)
c    V4.10.6   8-Jun-2005   AK' Inserted additional column 'ocltn' in .ret file giving occultation number
c    V4.10.6   9-Jun-2005   AK' Introduced a dynamic expression for all lats for smoothing
c    V4.10.7   6-Jul-2005   AK' Added feature to use vaa=vbar in case of CH3CN retrieval
c    V4.10.7   7-Jul-2005   AK' Added a priori constraint above the balloon for CH3CN retrieval
c    V4.11.0  16-Jan-2006   Adapted for new .ray and .mav output files.
c    V4.12.0  31-Jan-2006   Adopted double-precision of all slant-columns
c                           to avoid G77 problems wihen slerr > 1.0E+37
c    V4.13.0  06-Mar-2006   Adapted to new ray file with extra variable (bend)
c    V4.14.0  16-Aug-2006   Removed self-test (unnecessary with air in vmr file).
c    V4.14.1  22-Jun-2007   Various mods to allow up to 175 atmospheric levels (ACE)
c    V4.14.2   7-Nov-2007   
c    V4.14.3  22-Jun-2008   Mods to support 30 km level spacing (Titan)
c    V4.14.4  31-Aug-2009   Support new format of .lav and .ray files & C*57 specname
c    V4.14.12 11-Sep-2012   Increased mspeci to 160

      implicit none
      include "../ggg_const_params.f"
      include "../ggg_int_params.f"

      integer*4 hilay,i,ieof,ilay,iobs,it,j,jlay,kcol,kgas,
     & ldel,lmax,lmin,lnbc,lo,lowlay,ltop,lun_dc,lun_vmr,lunr_ray,
     & lunr_lav,lunw_ret,lunw_con,lunr_mav,lunr_rpt,lr,lunw_jac,nrow,ls,
     & lunr_vmr,lunw_vmr,nlhead,ncol_vmr,nlev_vmr,lvf,
     & malt,mcomm,mwin,mm,mobs,mode,msza,nalt,ilev,jgas,
     & mauxlav,nauxlav,ndrec,lrt,kk,
     & nauxmav,nauxret,ncollav,ncol_mav,ncomm,ngas,nl,nlev_ray,nlev_mav,
     & nn,nobs,noly,np,nsza,nspeci,ocltn,neq,ext
      parameter (lunr_mav=12)  ! runlog.mav file
      parameter (lun_vmr=13)  ! runlog.vmr output file (dormant) 
      parameter (lunr_ray=14)  ! runlog.ray  File of Slant Paths
      parameter (lunr_lav=15)  ! runlog.lav File of slant columns
      parameter (lunw_ret=16)  ! runlog.lav.ret  output file 
      parameter (lunw_con=21)  ! runlog.lav.ret  output file 
      parameter (lunr_vmr=22)  ! runlog.lav.ret  output file 
      parameter (lunw_vmr=23)  ! runlog.lav.ret  output file 
      parameter (lunr_rpt=17)  ! dirunret.rpt file 
      parameter (lun_dc=18)   ! diurnal.xxx input files 
      parameter (lunw_jac=20)  ! Jaconians for computing averaging kernels
      parameter (malt=40)     ! Max Number of levels in diurnal correction file
      parameter (mcomm=8)     ! Max Number of comment lines
      parameter (mwin=444)    ! Max Number of gases/windows to be retrieved
c      parameter (mlev=175)    ! Max Number of model layers
      parameter (mobs=310)    ! Max Number observations/spectra per occultation
      parameter (msza=130)    ! Max Number of angles in diurnal correction file
      parameter (mm=mobs+mlev)
      parameter (mauxlav=24)  ! Number of auxiliary columns in lav file
      parameter (nauxmav=4)   ! Number of auxiliary columns in mav file
      parameter (nauxret=12)  ! Number of auxiliary columns in ret file
      integer*4 jpvt1(mlev),jpvt2(mlev)
c
      character
     & col1*1,comment(mcomm-1)*1000,dcfile*(mfilepath),
     & dcflabel*1024,dcflag*2,
     & dcfstarr(msza)*8,gas1*32,lavlabel*12000,mavlabel*2000,
     & lavstarr(2*mwin+mauxlav)*32,ray_string*2400,specname*(nchar),
     & mavstarr(mspeci+nauxmav)*32,lavfile*40,spectrum_ray*(nchar),
     & specraywas*(nchar),start_spectrum_ray*(nchar),version*48,
     & next_spectrum_mav*(nchar),ray_format*20,
     & vmrfile*(mfilepath),modfile*80,gasnames(mgas+1)*32,c1*1,
     & gggdir*(mpath),dl*1,     !ggg directory path (GGGPATH?)
     & input_fmt*40,output_fmt*40
c
      logical dcfilexst
c
      real*8 aszai,aszaz,ddum,denom,doty(mobs),dotyi,dotyz,dsdot,dz,
     & delz,eorv,ervc,hour(mobs),houri,hourz,oblat(mobs),oblon(mobs),
     & solar_time,tlat(mobs),tlati,tlatz,tlon(mobs),tloni,tlonz,totwt,
     & totwz,totwz2,wt,year(mobs),yeari,yearz,zobs(mobs),mtlat,smfact
c
      real*4 afit,pins(mobs),tins(mobs),asza_lav(mobs),azim(mobs),
     & cfit,d(mlev),d2r,dcf,
     & vmr(mlev,mgas),fuap,rr2,
     & dnight,dsp(mobs),dsunny,
     & dum,dcell,erad,ev(mlev,mwin),fr,frac,asza_ray,
     & p(mlev),pd1(mm,mlev+1),pd2(mm,mlev+1),pobs(mobs),
     & rnorm,run(mobs),scht,slsfin, smoo,sp(mobs,mlev),sza,t(mlev),ta,
     & tab(malt,msza),talt(malt),tc,tinfo,tobs(mobs),hobs(mobs),
     & tsza(msza),tt,ttv,ttx,tu,tv,twodint,tx,ufit,vaa(mlev),vbar,
     & vcon(mlev,mwin),dzmin,dzmax,zminwas,
     & vold,vunc(mlev,mwin),wk(mlev),xx,z(mlev),zeff(mobs),
     & zmin(mobs),zz(mm),verr,tverr,km2cm

      real*8
     &  slcol(mobs,mwin),slerr(mobs,mwin),vaacol,concol,unccol

      parameter(d2r=spi/180,erad=6378.,scht=5.5, km2cm=100000.)
c
c     Platform specification:      DG090519
      call get_ggg_environment(gggdir, dl)
      lrt=lnbc(gggdir)     !Length of gggdir
c

      ray_format='(a,7f11.5,200f11.5)'
      version=' DIURNRET     Version 4.16      2013-01-31   GCT'
      open(lunr_rpt,file='diurnret.rpt',status='unknown')
      write(6,*)version
      write(lunr_rpt,*)version
c
      zminwas=0.0  ! avoid compiler warnings
      dzmax=0.0  ! avoid compiler warnings
      if (iargc() == 0 ) then
         write(6,117)
 117     format(' Enter Slant Column File (e.g. fai97rat.lav): ',$)
         read(5,'(a)')lavfile
      elseif (iargc() == 1) then
         call getarg(1, lavfile)
      else
         stop 'Usage: $gggpath/bin/diurnret lavfile'
      endif
      write(lunr_rpt,*) lavfile
      lo=lnbc(lavfile)
c==================================================================
      write(*,*)
      ieof=0
c  Read slant paths, slant columns, and slant column errors for each spectrum
c      write(*,*)lo,lavfile(:lo-4)//'.mav'
      open(lunr_mav,file=lavfile(:lo-4)//'.mav',status='old')
      read(lunr_mav,*)   ! read version string
      open(lunr_ray,file=lavfile(:lo-4)//'.ray',status='old')
      read(lunr_ray,*)nn,nlev_ray
      nlev_ray=nlev_ray-7-ncell  ! First 7 cols are Spectrum Zobs Pres SZA Bend FOV Zmin
      if(nlev_ray.gt.mlev) stop 'Increase parameter MLEV'
      do i=2,nn+1
        read(lunr_ray,'(a)') ray_string
        if(ray_string(:7).eq.'format=') ray_format=ray_string(8:)
      end do
      open(lunr_lav,file=lavfile,status='old')
      read(lunr_lav,'(i2,i4,i7,i4)')ncomm,ncollav,nrow,nauxlav

      if(nauxlav.gt.mauxlav) stop 'nauxlav .gt. mauxlav'
      if(ncomm.gt.mcomm) write(6,*)'Warning: Increase parameter MCOMM'
      do i=2,ncomm-1
         read(lunr_lav,'(a)') comment(i)
      end do
      read(lunr_lav,'(a)')lavlabel
      call substr(lavlabel,lavstarr,2*mwin+nauxlav,nn)
      if(nn.ne.ncollav) then
         write(*,*)'# columns according to file header: ',ncollav
         write(*,*)'# of column labels actually found: ',nn
      stop 'Inconsistent number of columns/labels'
      endif
      ngas=(ncollav-nauxlav)/2
      if(ngas.gt.mwin) then
         write(6,*) 'Increase parameter MGAS to',ngas
         stop
      endif
      if( index(lavlabel,'spectrum').gt.0 ) then
         ls=57
         input_fmt='(a57,f13.8,NNf13.5,MMM(e12.4))'
         write(input_fmt(12:13),'(i2.2)') nauxlav-2
         write(input_fmt(20:22),'(i3.3)') ncollav-nauxlav
      else
         ls=1
         input_fmt='(a1,f13.8,NNf13.5,MMM(e12.4))'
         write(input_fmt(11:12),'(i2.2)') nauxlav-1
         write(input_fmt(19:21),'(i3.3)') ncollav-nauxlav
      endif
c     write(*,*) 'input format =  ',input_fmt
c      input_fmt = comment(ncomm-1)(:40)
c      write(*,*) 'input format =  ',input_fmt

c Write .jac file
      open(lunw_jac,file=lavfile(:lo)//'.jac',status='unknown')
c
c Write .ret file
      ndrec=0
      open(lunw_ret,file=lavfile(:lo)//'.ret',status='unknown')
      open(lunw_con,file=lavfile(:lo)//'.con',status='unknown')
      write(lunw_ret,*)ncomm+1,nauxret+2*ngas  
      write(lunw_con,*)ncomm+1,nauxret+2*ngas  
      write(lunw_ret,*)version
      write(lunw_con,*)version
      do i=2,ncomm-2
         write(lunw_ret,'(a)') comment(i)(:lnbc(comment(i)))
         write(lunw_con,'(a)') comment(i)(:lnbc(comment(i)))
      end do
      output_fmt='(i6,f12.6,7f12.5,3(1pe12.4),MMMe12.4)'
      write(output_fmt(29:31),'(i3.3)') nauxret+2*ngas
      write(lunw_ret,'(a)') output_fmt
      write(lunw_con,'(a)') output_fmt
      gas1=lavstarr(nauxlav+1)
      write(lunw_ret,'(a)')
     &' ocltn  year          day           hour       tplat'//
     &'      tplong      asza         altitude'//
     &'   temp       pres        dens       theta      '//
     & lavlabel(index(lavlabel,gas1(:lnbc(gas1)))-2:lnbc(lavlabel)+1)
      write(lunw_con,'(a)')
     &' ocltn  year          day           hour       tplat'//
     &'      tplong      asza         altitude'//
     &'   temp       pres        dens       theta      '//
     & lavlabel(index(lavlabel,gas1(:lnbc(gas1)))-2:lnbc(lavlabel)+1)
c
c  Main loop over occultations
      do ocltn=1,999999
      lowlay=nlev_ray  ! lowest atmospheric layer used
      hilay=0          ! highest unused atmospheric layer
      do  iobs=1,mobs  ! Loop over good spectra within one occultation
 2      continue
c        if( index(lavlabel,'spectrum').gt.0 ) then
c        if(spectrum_flag.eq.1) then
           read(lunr_lav,input_fmt,end=14)
     &     specname(:ls),year(iobs),doty(iobs),hour(iobs),
     &     run(iobs),oblat(iobs),oblon(iobs),zobs(iobs),zmin(iobs),
     &     asza_lav(iobs),azim(iobs),
     &     dum,dum,dum,dum,tins(iobs),pins(iobs),
     &     tobs(iobs),pobs(iobs),hobs(iobs),dum,dum,dum,dum,
     &     (slcol(iobs,kgas),slerr(iobs,kgas),kgas=1,ngas)
c        else
c           read(lunr_lav,input_fmt,end=14)
c     &     specname(:ls),year(iobs),doty(iobs),hour(iobs),
c     &     run(iobs),oblat(iobs),oblon(iobs),zobs(iobs),zmin(iobs),
c     &     asza(iobs),azim(iobs),dum,dum,dum,dum,tins(iobs),pins(iobs),
c     &     tobs(iobs),pobs(iobs),hobs(iobs),dum,dum,dum,dum,
c     &     (slcol(iobs,kgas),slerr(iobs,kgas),kgas=1,ngas)
c        endif
        col1=specname(ls:ls)
c
c  Find largest tangent altitude separation
        if(iobs.gt.1) then
           dzmin=abs(zmin(iobs)-zminwas)
           if(dzmin.gt.dzmax) dzmax=dzmin
        endif
        zminwas=zmin(iobs)
c
c        write(*,*) ray_format
c        write(*,*) ray_string
        read(ray_string,ray_format) spectrum_ray,dum,dum,
     &  asza_ray,dum,dum,dum,(dcell,j=1,ncell),
     &  (sp(iobs,jlay),jlay=1,nlev_ray)
        read(lunr_ray,'(a)',end=16) ray_string
        lr=lnbc(spectrum_ray)
c        write(*,'(a)')col1//spectrum_ray(:lr+3)//specraywas(:lr+3)//next_spectrum_mav
c     &//ray_string(:lr)
        if (ray_string(lr-3:lr).eq.'.adf') then
c          write(*,*)'ACE spectra'
           ext=7 ! ACE spectra all have .adf ending, so you need to match .XXX.adf
        else
           ext=3
        endif
        do while(ray_string(lr-ext:lr).eq.spectrum_ray(lr-ext:lr))
           read(ray_string,ray_format) spectrum_ray,dum,dum,
     &     asza_ray,dum,dum,dum,(dcell,j=1,ncell),
     &     (sp(iobs,jlay),jlay=1,nlev_ray)
           read(lunr_ray,'(a)',end=16) ray_string
        end do
16      if(col1.eq.'-' .or. col1.eq.':') go to 2  ! skip bad spectra
c
        if(asza_lav(iobs).ne.asza_ray) write(6,'(a,i4,2f10.4)')
     &   spectrum_ray(:24)//ray_string(:24)//
     &  ': Warning: zenith angle mismatch ',iobs,asza_lav(iobs),asza_ray
c
c  Compute Latitude and longitude of tangent points
        if(asza_lav(iobs).le.90.) then
           zeff(iobs)=zobs(iobs)+scht*cos(d2r*asza_lav(iobs))
           tlat(iobs)=oblat(iobs)
           tlon(iobs)=oblon(iobs)
        else
           zeff(iobs)=zmin(iobs)
           call zenaz(2,oblat(iobs),oblon(iobs),zobs(iobs),
     &     int(year(iobs)),1,int(doty(iobs)),hour(iobs)/24,
     &     ddum,ddum,eorv,ervc,tlat(iobs),tlon(iobs),ddum)
        endif
c
        do ilay=1,nlev_ray
           if(sp(iobs,ilay).gt.0.001) go to 76
        end do
        write(*,*) 'Warning: slant paths are all zero for ',spectrum_ray
 76     if(ilay.lt.lowlay) lowlay=ilay  ! lowest atmospheric layer used
        if(ilay.gt. hilay) hilay =ilay  ! highest unused atmospheric layer
 
        if(iobs.le.1 ) start_spectrum_ray=spectrum_ray
        specraywas=spectrum_ray
        if(col1.eq.'+') go to 15
      end do   !  iobs=1,mobs
      read(lunr_lav,*,end=14)
      stop ' Increase parameter MOBS'
 14   ieof=1      ! EOF
      iobs=iobs-1
 15   nobs=iobs
c      write(*,*)'ocltn,nobs,ieof= ',ocltn,nobs,ieof
      if(nobs .le. 0) then
         if(ieof.gt.0 ) then
            go to 99
         else
            stop 'nobs=0'
         endif
      endif
      nl=nlev_ray-lowlay+1
c      write(*,*)'nlev_ray, lowlay, nl ',nlev_ray, lowlay, nl
c      write(*,*)'Next Spectrum=',next_spectrum_mav,spectrum_ray
c==================================================================
c  Read model (altitudes, temperatures, pressures, & densities).
c  The vmrs, which are in the same file, are not used in the retrieval.
c
c Only do this if the day has changed. If we are retrieving profiles
c from ascent & sunset spectra acquired on the same day, don't read mavfile.
c      write(*,'(a)') spectrum_ray(4:8)//'  '//next_spectrum_mav(4:8)
      do while( spectrum_ray(4:8) .ne. next_spectrum_mav(4:8) )
         write(6,*)
         read(lunr_mav,'(14x,a)',end=99) next_spectrum_mav 
         write(*,*) 'Next Spectrum =',next_spectrum_mav
         read(lunr_mav,*)nn,ncol_mav,nlev_mav
         nlev_mav=nlev_mav-ncell    ! ignore cell levels
         if(nlev_mav.ne.nlev_ray) stop 'nlev_mav .ne. nlev_ray'
         nspeci=ncol_mav-nauxmav   ! First 4 columns are Z T P D
         if(nspeci.gt.mspeci) then
            write(6,*)' mspeci, nspeci = ',mspeci,nspeci
            write(6,*)' ncol_mav, nauxmav = ',ncol_mav,nauxmav
            stop 'Increase parameter MSPECI'
         endif
         do i=1,nn-4
           read(lunr_mav,*)  ! skip unecessary indo (e.g. trop altitude)
         end do
         read(lunr_mav,'(a)') vmrfile
         read(lunr_mav,'(a)') modfile
         read(lunr_mav,'(a)') mavlabel  ! column labels
         call lowercase(mavlabel)
         call substr(mavlabel,mavstarr,mspeci+nauxmav,kcol)
         if( kcol .ne. ncol_mav ) then
            write(*,*)mavlabel(:lnbc(mavlabel))
            write(*,*)'kcol,ncol_mav=',kcol,ncol_mav
            stop 'column # mismatch in mavfile'
         endif
c
c  Skip cell info
         do i=1,ncell
           read(lunr_mav,'(2f7.2,2e11.3)')z(1),t(1),p(1),d(1)
         end do
c  Read atmospheric level info
         do jlay=1,nlev_mav
          read(lunr_mav,'(2f7.2,2e11.3)')z(jlay),t(jlay),p(jlay),d(jlay)
         end do
      end do      !   spectrum_ray(4:8) .ne. next_spectrum_mav(4:8) 
c      write(6,'(i4,2x,a)') ocltn,specraywas(4:9)
      write(lunr_rpt,*) specraywas(4:9)
      write(6,*)       ' Number of different retrieved gases= ',ngas
      write(lunr_rpt,*)' Number of different retrieved gases= ',ngas
      write(6,*)       ' Number of observed spectra used    = ',nobs
      write(lunr_rpt,*)' Number of observed spectra used    = ',nobs
      write(6,*)       ' Lowest atmospheric layer used      = ',lowlay
      write(lunr_rpt,*)' Lowest atmospheric layer used      = ',lowlay
      write(6,*)       ' Highest atmospheric tangent layer  = ',hilay+1
      write(lunr_rpt,*)' Highest atmospheric tangent layer  = ',hilay+1
      write(6,*)       ' Number of retrieved layers         = ',nlev_ray
      write(lunr_rpt,*)' Number of retrieved layers         = ',nlev_ray
      write(lunr_rpt,*)' GAS       Vbar    Vmda    A_Fit(%) '
     &//' C_Fit(%) U_Fit(%)  Frac(%)  SMOO'
c==================================================================
c  Read VMR file
      lvf=lnbc(vmrfile)
      write(*,*) 'Reading ',vmrfile
      open(lunr_vmr,file=vmrfile,status='old')
      read(lunr_vmr,*) nlhead,ncol_vmr
      do j=1,nlhead-1
         read(lunr_vmr,'(a)') comment(j)
      end do
      call substr(comment(nlhead-1), gasnames, mgas,kk)
      if(kk.ne.ncol_vmr) stop 'kk.ne.ncol_vmr'
      do ilev=1,mlev
      read(lunr_vmr,*,end=78) z(ilev),(vmr(ilev,jgas),jgas=1,ncol_vmr-1)
      end do
78    nlev_vmr=ilev-1
      close(lunr_vmr)

c-------------------------------------------------------------------
c  Start retrievals, gas by gas
      do 300 kgas=1,ngas  ! Main Loop over retrieval gases
c  Set up the equation  YY = PD.X  where YY=slcol/slerr & PD=dens*slpath/slerr
c  Also compute mean vmr (VBAR) for profile.
        do jgas=1,ncol_vmr-1
           c1=lavstarr(nauxlav+2*kgas-1)(1:1)
           if(c1.eq.'1' .or. c1.eq.'t') then
               lavstarr(nauxlav+2*kgas-1)=lavstarr(nauxlav+2*kgas-1)(2:)
           endif
           call lowercase(gasnames(jgas+1))
           if(lavstarr(nauxlav+2*kgas-1).eq.gasnames(1+jgas)) exit
        end do
        write(*,*) kgas,jgas,lavstarr(nauxlav+2*kgas-1),gasnames(jgas+1)
        
        tinfo=0.0
        do iobs=1,nobs
          dsp(iobs)=0.0
          do j=1,nl
             xx=d(j+lowlay-1)*km2cm*sp(iobs,j+lowlay-1)
             dsp(iobs)=dsp(iobs)+xx
c              write(*,*)'xx slerr=',xx,slerr(iobs,kgas)
             pd1(iobs,j)=xx/slerr(iobs,kgas)
             if(slerr(iobs,kgas).gt.1.E+36)
     &       write(*,*)' slerr.gt.1.0e+36',
     &       iobs,kgas,j,xx,slerr(iobs,kgas),pd1(iobs,j),
     &       specname,lavstarr(nauxlav+2*kgas-1)
          end do
          pd1(iobs,nl+1)=slcol(iobs,kgas)/slerr(iobs,kgas)
          tinfo=tinfo+pd1(iobs,nl+1)**2
        end do ! iobs=1,nobs
        tinfo=sqrt(tinfo)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  Multiply array PD1 by diurnal correction factors, if appropriate
        solar_time=(hour(1)+tlon(1)/15 + hour(nobs)+tlon(nobs)/15)/2
        if( dmod(solar_time+24,24.d0) .lt. 9.0) then
           dcflag='sr'
        elseif( dmod(solar_time+24,24.d0) .gt. 13.5) then
           dcflag='ss'
        else
           dcflag='  ' ! don't diurnally correct data acquired around noon
        endif
        dcfile=gggdir(:lrt)//'diurnal/'//vmrfile(lvf-8:lvf-4)//'/dcf.'//
     & dcflag//'.'//lavstarr(2*kgas+nauxlav-1)
        inquire(file=dcfile,exist=dcfilexst)
        if(dcfilexst) then
           write(*,*) 'Found dcfile=',dcfile
           open(lun_dc,file=dcfile,status='old')
           read(lun_dc,*) nn, nsza
           do i=2,nn
              read(lun_dc,'(a)')dcflabel
           end do
c           call skiprec(lun_dc,nn-2)
c           read(lun_dc,'(a)')dcflabel
           call substr(dcflabel,dcfstarr,msza,nn)
           if(nn .ne. nsza) then
              write(*,*) nn, nsza
              stop 'column/label mismatch'
           endif
           nsza=nsza-1    ! first column is the altitudes
           read(dcflabel(9:),*)(tsza(j),j=1,nsza)  ! solar zenith angles
           do nalt=1,malt
              read(lun_dc,*,end=11)talt(nalt),(tab(nalt,j),j=1,nsza)
           end do
           stop 'increase parameter MALT'
 11        close(lun_dc)
           nalt=nalt-1
c
c  Correct any -ve solar zenith angles (used to denote am in RJS model)
           do j=1,nsza
              tsza(j)=abs(tsza(j))
           end do
c
c Average DCF's from the near (night) and far (sunny) sides of tangent point
           do iobs=1,nobs
              do j=1,nl
                 if(z(j).gt.zmin(iobs)) then
                    if(zmin(iobs).eq.zobs(iobs)) then  !  up looking geometry
                     sza=asin((zobs(iobs)+erad)*sin(d2r*asza_lav(iobs))/
     &                 (z(j)+erad))/d2r
                    else                   !  limb viewing geometry
                       sza=asin((zmin(iobs)+erad)/(z(j)+erad))/d2r
                    endif
                   dsunny=twodint(z(j),sza,talt,tsza,tab,malt,nalt,nsza)
                    dnight=twodint(z(j),180-sza,talt,tsza,tab,malt,nalt,
     &              nsza)
                    fr=(zobs(iobs)-z(j))/(z(j+1)-z(j))
                    if(fr.le.-0.5) then     ! levels well above observer
                       dcf=dsunny
                    elseif(fr.ge.+0.5) then ! levels well below observer e.g. ATMOS
                       dcf=(dsunny+dnight)/2
                    else                    ! levels immediately  adjacent to ZOBS
                       dcf=((3-2*fr)*dsunny+(1+2*fr)*dnight)/4
                    endif
                    pd1(iobs,j)=pd1(iobs,j)*dcf
                 endif   ! (z(j).gt.zmin(iobs))
              end do  ! j=1,nl
           end do  ! iobs=1,nobs
        else
           dcflag='  '   !  diurnal correction file not found
           write(*,*) 'Not found dcfile=',dcfile
        endif   ! (dcfilexst)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  Compute an approximate vmr profile (VAA) for the purpose of assigning
c  relative weights to the constraints.
        ttv=0.0
        ttx=0.0
        do j=1,nl
          tv=0.0
          tx=0.0
          do iobs=1,nobs
            wt=abs(d(j+lowlay-1)*km2cm*sp(iobs,j+lowlay-1)/
     &      slerr(iobs,kgas))
            tv=tv+wt*(slcol(iobs,kgas)/dsp(iobs))
            tx=tx+wt
          end do
          vaa(j)=tv/tx
          ttv=ttv+tv
          ttx=ttx+tx
        end do
        vbar=abs(ttv/ttx)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c Smooth VAA twice (this should not be necessary but it is)
        do it=1,2
           vold=vaa(1)
           do j=2,nl-1
              vaa(j)=(vold+vaa(j)+vaa(j+1) )/3
              vold=vaa(j)
           end do
        end do
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  Scale PD's by vbar to avoid overflows. 
        do iobs=1,nobs
           do j=1,nl
              pd1(iobs,j)=pd1(iobs,j)*(abs(vaa(j))+vbar)/2
           end do
        end do
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  Augment matrix of partials (PD) with NLAY derivative constraints (smoothing)
c  Augment measurement vector (YY) with NLAY zeros
        do ilay=1,nl
           jpvt1(ilay)=0
           jpvt2(ilay)=0
           do jlay=1,nl
              pd1(nobs+ilay,jlay)=0.0
           end do
           pd1(nobs+ilay,nl+1)=0.0
        end do
C
C  Choose the smoothing coefficient SMOO so that it always supplies a
c  significantly lesser constraint than the measurements themselves.
c  This way, it does not significantly impact (reduce) the error bars.
         mtlat=(tlat(1)+tlat(nobs))*dpi/360 ! in rad
c         smfact=0.015  ! MATMOS
c         smfact=0.0625*cos(mtlat)+0.01*(1-cos(mtlat))**2  ! MkIV balloon
         smfact=0.05*cos(mtlat)+0.01*(1-cos(mtlat))**2  ! MkIV balloon
         smfact=0.020  !  for ACE tropical data
         smfact=0.010  ! VATMOS
         smoo=sqrt(smfact*tinfo) ! dynamic expression for all latitudes
c        smoo=sqrt(tinfo/16)  ! Normal (mid-latitudes)
c        smoo=sqrt(tinfo/64)  ! for N2O isotopes
c        smoo=sqrt(tinfo/32)  ! for Esrange profiles
c        smoo=sqrt(tinfo/8)   ! for HOCl and H2O2
c        smoo=6.              ! Old way (prior to V320)

c original smoothing criterion
        neq=nobs+nl ! no. of equations
        pd1(nobs+1,1)=smoo
        pd1(nobs+1,2)=-smoo
        do ilay=2,nl-1
           noly=nobs+ilay
           pd1(noly,ilay-1)=-smoo
           pd1(noly,ilay)=2*smoo
           pd1(noly,ilay+1)=-smoo
        end do
        pd1(nobs+nl,nl)=smoo
        pd1(nobs+nl,nl-1)=-smoo

c   Duplicate entire array of partial differentials
        write(lunw_jac,'(4i5)')ocltn,kgas,nobs,nl
        do iobs=1,neq
           do ilay=1,nl+1
              pd2(iobs,ilay)=pd1(iobs,ilay)
           end do
           write(lunw_jac,'(200(1pe12.4))')(pd1(iobs,jlay),jlay=1,nl+1)
        end do
c
c  Perform actual retrievals; the first constrained, the second unconstrained
        call nnls(pd1,mm,neq,nl,pd1(1,nl+1),vcon(1,kgas),
     &  rnorm,wk,zz,jpvt1,mode)
        call snnls(pd2,mm,neq,nl,vunc(1,kgas),np,jpvt2,wk,0)
        rnorm=slsfin(pd2,mm,neq,nl,np,jpvt2)
c        write(*,*) 'lowlay=',lowlay,vmr(3,jgas),vmr(4,jgas),vmr(5,jgas)
        fuap=0.25  ! Fractional uncertainty in A Priori vmr profiles
        do ilay=1,nl
          xx=(abs(vaa(ilay))+vbar)/2
          vcon(ilay,kgas)=xx*vcon(ilay,kgas)
          vunc(ilay,kgas)=xx*vunc(ilay,kgas)
          ev(ilay,kgas)=
     &    xx*dsqrt(dsdot(nl,pd2(ilay,1),mm,pd2(ilay,1),mm))
          if(jgas.ne.2 .and. jgas.ne.41) then
          rr2=(ev(ilay,kgas)/(fuap*vmr(lowlay+ilay-1,jgas)))**2
          vmr(lowlay+ilay-1,jgas)=
     & (vcon(ilay,kgas)+ vmr(lowlay+ilay-1,jgas)*rr2)/(1+rr2)
c     & (vcon(ilay,kgas)/ev(ilay,kgas)**2+
c     & 1/vmr(lowlay+ilay-1,jgas)*fuap**2)/
c     &   (1/ev(ilay,kgas)**2 + 1/(vmr(lowlay+ilay-1,jgas)*fuap)**2)
          end if
        end do

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  Recalculate slant columns using original and both constrained (VCON) &
c  unconstrained (VUN)  retrieved vmr profiles.
c  Compute the "goodness of fit" to the measured slant columns.
c  Compute the factor (FRAC) by which the deviations of the measurements
c  from the re-calculated slant columns exceeds their original error estimates.
        tverr=0.0
        ta=0.0
        tc=0.0
        tu=0.0
        tt=0.0
        do iobs=1,nobs
          vaacol=0.0
          concol=0.0
          unccol=0.0
          do ilay=1,nl
            xx=d(ilay+lowlay-1)*km2cm*sp(iobs,ilay+lowlay-1)
            vaacol=vaacol+xx*vaa(ilay)
            unccol=unccol+xx*vunc(ilay,kgas)
            concol=concol+xx*vcon(ilay,kgas)
          end do
c          write(6,*)slcol(iobs,kgas),concol,unccol
          tverr=tverr+(vaacol/slerr(iobs,kgas))**2
          ta=ta+((slcol(iobs,kgas)-vaacol)/slerr(iobs,kgas))**2
          tc=tc+((slcol(iobs,kgas)-concol)/slerr(iobs,kgas))**2
          tu=tu+((slcol(iobs,kgas)-unccol)/slerr(iobs,kgas))**2
          tt=tt+ (slcol(iobs,kgas)/slerr(iobs,kgas))**2
c          write(6,*)ta,tc,tu,tt
        end do
        afit=sqrt(ta/tt)
        cfit=sqrt(tc/tt)
        ufit=sqrt(tu/tt)
        frac=sqrt(tu)/nobs
        verr=vbar/sqrt(tverr)
      write(lunr_rpt,'(a10,1p2e9.2,0p4f8.2,0pf7.1,a2)')
     &  lavstarr(2*kgas+nauxlav-1),vbar,
     &  verr,100*afit,100.*cfit,100*ufit,100*frac,smoo,dcflag
 300  continue    ! kgas=1,ngas   Main Loop over retrieval gases
c==================================================================
c  Write new VMR file
      lvf=lnbc(vmrfile)
      write(*,*) 'Reading ',vmrfile
      open(lunw_vmr,file=vmrfile(:lvf)//'.new',status='unknown')
      write(lunw_vmr,*) nlhead,ncol_vmr
      do j=1,nlhead-1
         write(lunw_vmr,'(a)') comment(j)(:lnbc(comment(j)))
      end do
      do ilev=1,nlev_vmr
         write(lunw_vmr,'(f7.1,200(1pe11.3))') z(ilev),
     &   (vmr(ilev,jgas),jgas=1,ncol_vmr-1)
      end do
      close(lunw_vmr)
c-----------------------------------------------------------------
c  Write out retrieved vmrs.
c  Create matrix of retrieved vmrs in a format which can be plotted by GPLOTEM
c  Write them in chronological order.
      ltop=hilay-lowlay+1+1
c      write(*,*)nobs,zmin(1),zmin(nobs),year(1),year(nobs)
      if( (year(nobs)-year(1))/(zmin(nobs)-zmin(1)) .gt. 0.0 ) then
        lmin=2
        lmax=ltop
        ldel=1
      else
        lmin=ltop
        lmax=2
        ldel=-1
      endif

      delz=(z(lmax)-z(lmin))/(lmax-lmin)  ! average level spacing
      do ilay=lmin,lmax,ldel
c
c  Interpolate YEAR, DAY, UT, ASZA, TLAT & TLON to the retrieval altitudes.
        totwt=0.0d0
        totwz=0.0d0
        totwz2=0.0d0
        yeari=0.0d0
        yearz=0.0d0
        dotyi=0.0d0
        dotyz=0.0d0
        houri=0.0d0
        hourz=0.0d0
        tlati=0.0d0
        tlatz=0.0d0
        tloni=0.0d0
        tlonz=0.0d0
        aszai=0.0d0
        aszaz=0.0d0
        do iobs=1,nobs
          dz=z(ilay+lowlay-1)-zeff(iobs)
c          write(*,*)iobs,year(iobs),z(ilay+lowlay-1),zeff(iobs),dz,delz
          wt=exp(-0.5*(dz/delz)**2)
          totwt=totwt+wt
          totwz=totwz+wt*dz
          totwz2=totwz2+wt*dz**2
          yeari=yeari+wt*year(iobs)
          yearz=yearz+wt*year(iobs)*dz
          dotyi=dotyi+wt*doty(iobs)
          dotyz=dotyz+wt*doty(iobs)*dz
          houri=houri+wt*hour(iobs)
          hourz=hourz+wt*hour(iobs)*dz
          tlati=tlati+wt*tlat(iobs)
          tlatz=tlatz+wt*tlat(iobs)*dz
          tloni=tloni+wt*tlon(iobs)
          tlonz=tlonz+wt*tlon(iobs)*dz
          aszai=aszai+wt*asza_lav(iobs)
          aszaz=aszaz+wt*asza_lav(iobs)*dz
        end do
        denom=(totwt*totwz2-totwz**2)
c        write(*,*)yeari,yearz,totwz2,totwz,totwt,denom
        yeari=(yeari*totwz2-yearz*totwz)/denom
        dotyi=(dotyi*totwz2-dotyz*totwz)/denom
        houri=(houri*totwz2-hourz*totwz)/denom
        tlati=(tlati*totwz2-tlatz*totwz)/denom
        tloni=(tloni*totwz2-tlonz*totwz)/denom
        aszai=(aszai*totwz2-aszaz*totwz)/denom
c
c        if(yeari.eq.2048.0) stop
        if(aszai.gt.90.) aszai=90.0

        ndrec=ndrec+1
        write(lunw_ret,output_fmt)
     &  ocltn,yeari,dotyi,houri,tlati,tloni,aszai,
     &  z(ilay+lowlay-1),t(ilay+lowlay-1),1013.25*p(ilay+lowlay-1),
     &  d(ilay+lowlay-1),
     &  t(ilay+lowlay-1)*p(ilay+lowlay-1)**(-.286),
     &  (vunc(ilay,kgas),ev(ilay,kgas),kgas=1,ngas)
        write(lunw_con,output_fmt)
     &  ocltn,yeari,dotyi,houri,tlati,tloni,aszai,
     &  z(ilay+lowlay-1),t(ilay+lowlay-1),1013.25*p(ilay+lowlay-1),
     &  d(ilay+lowlay-1),
     &  t(ilay+lowlay-1)*p(ilay+lowlay-1)**(-.286),
     &  (vcon(ilay,kgas),ev(ilay,kgas),kgas=1,ngas)
c
c  PD1 contains the constrained retrieved vmrs 
c  PD2 contains the unconstrained retrieved vmrs
c  SP  contains their uncertainties.
      end do   ! ilay=2,hilay-lowlay+1+3
c=======================================================================
      if(ieof.gt.0) go to 99
      end do  ! ocltn=1,999999
99    close(lunr_ray)
      close(lunr_mav)
      close(lunr_lav)
      close(lunr_rpt)
      close(lunw_jac)
      close(lunw_ret)
      close(lunw_con)
      end
