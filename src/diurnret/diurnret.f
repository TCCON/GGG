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
c  ' V3.0.1 23-Mar-96   BS' increased mvmr from 62 to 65
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
c  ' V4.9.4  7-Sep-99  GCT' Increased mvmr from 65 to 122 for isotopic retieval
c  ' V4.9.5  8-Sep-99  GCT' Uses column labels from .mav file (instead of MOLNUM) to identify gases.
c  ' V4.9.7  5-Jan-00  GCT' Increased mvmr from 100 to 122 for isotopic retieval
c  ' V4.9.8 28-May-00  GCT' Defined LTOP variable to simplify code.
c  ' V4.9.9 19-Jul-00  GCT' Added IF's to open right input files for ATMOS
c  ' V4.9.9  9-Jan-01  GCT' Increased MCOMM to 6
c    V4.9.10  29-Oct-2001   GCT' Added "d0" to DMOD statements. 9000 format
c    V4.10.1  29-Jul-2002   Calculate Minimum Detectable Abundance vmr (Vmda)
c    V4.10.2  30-Jul-2002   GCT' Added "diurnret.out" output file.
c    V4.10.3  23-Jul-2003   GCT' Increased parameter msza from 22 to 130
c    V4.10.3  23-Jul-2003   Increased dimension of dcflabel from 256 to 1024
c    V4.10.4   1-Nov-2004   Increased parameter mvmr from 126 to 140
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

      implicit none
      include "../ggg_const_params.f"
      include "../ggg_int_params.f"

      integer*4 hilay,i,jcol,ieof,ilay,iobs,it,ix,ig,j,jlay,kcol,kgas,
     & ldel,lmax,lmin,lnbc,lo,lowlay,ltop,lun_dc,lun_vmr,lun_ray,
     & lun_lav,lun_ret,lun_mav,lun_rpt,lr,lun_jac,nrow,ls,
     & malt,mcomm,mgas,mlay,mm,mobs,mode,molid,msza,mvmr,nalt,
     & mauxlav,nauxlav,ndrec,
     & nauxmav,nauxret,ncollav,ncol_mav,ncomm,ngas,nl,nlev_ray,nlev_mav,
     & nn,nobs,noly,np,nsza,nvmr,ocltn,nlab,neq,ext
      parameter (lun_mav=12)  ! runlog.mav file
      parameter (lun_vmr=13)  ! runlog.vmr output file (dormant) 
      parameter (lun_ray=14)  ! runlog.ray  File of Slant Paths
      parameter (lun_lav=15)  ! runlog.lav File of slant columns
      parameter (lun_ret=16)  ! runlog.lav.ret  output file 
      parameter (lun_rpt=17)  ! dirunret.rpt file 
      parameter (lun_dc=18)   ! diurnal.xxx input files 
      parameter (lun_jac=20)  ! Jaconians for computing averaging kernels
      parameter (malt=40)     ! Max Number of levels in diurnal correction file
      parameter (mcomm=8)     ! Max Number of comment lines
      parameter (mgas=444)    ! Max Number of gases/windows to be retrieved
      parameter (mlay=175)    ! Max Number of model layers
      parameter (mobs=310)    ! Max Number observations/spectra per occultation
      parameter (msza=130)    ! Max Number of angles in diurnal correction file
      parameter (mvmr=150)    ! Max Number of vmr profiles to be handled
      parameter (mm=mobs+mlay)
      parameter (mauxlav=24)  ! Number of auxiliary columns in lav file
      parameter (nauxmav=4)   ! Number of auxiliary columns in mav file
      parameter (nauxret=12)  ! Number of auxiliary columns in ret file
      integer*4 jpvt1(mlay),jpvt2(mlay)
c
      character
     & col1*1,comment(mcomm-1)*80,dcfile*48,dcflabel*1024,dcflag*2,
     & dcfstarr(msza)*8,gas*10,gas1*32,lavlabel*12000,mavlabel*1680,
     & lavstarr(2*mgas+mauxlav)*32,ray_string*2400,specname*57,
     & mavstarr(mvmr+nauxmav)*32,lavfile*40,runlab*57,runwas*57,
     & start_runlab*57,version*44,next_spectrum*34,ray_format*20,
     & input_fmt*40,output_fmt*40
c
      logical dcfilexst
c
      real*8 aszai,aszaz,ddum,denom,doty(mobs),dotyi,dotyz,dsdot,dz,
     & delz,eorv,ervc,hour(mobs),houri,hourz,oblat(mobs),oblon(mobs),
     & solar_time,tlat(mobs),tlati,tlatz,tlon(mobs),tloni,tlonz,totwt,
     & totwz,totwz2,wt,year(mobs),yeari,yearz,zobs(mobs),mtlat,smfact,
     & aprierr
c
      real*4 afit,pins(mobs),tins(mobs),asza(mobs),azim(mobs),
     & cfit,d(mlay),d2r,dcf,
     & dnight,dsp(mobs),dsunny,
     & dum,dcell,erad,ev(mlay,mgas),fr,frac,obs_asza,
     & ofit,p(mlay),pd1(mm,mlay+1),pd2(mm,mlay+1),pobs(mobs),
     & rnorm,run(mobs),scht,slsfin, smoo,sp(mobs,mlay),sza,t(mlay),ta,
     & tab(malt,msza),talt(malt),tc,tinfo,to,tobs(mobs),hobs(mobs),
     & tsza(msza),tt,ttv,ttx,tu,tv,twodint,tx,ufit,vaa(mlay),vbar,
     & vcon(mlay,mgas),vmin,vmr(mvmr,mlay),dzmin,dzmax,zminwas,
     & vold,vunc(mlay,mgas),wk(mlay),xx,z(mlay),zeff(mobs),
     & zmin(mobs),zz(mm),verr,tverr,km2cm

      real*8
     &  slcol(mobs,mgas),slerr(mobs,mgas),vaacol,orgcol,concol,unccol

      parameter(d2r=spi/180,erad=6378.,scht=5.5, km2cm=100000.)
c
      ray_format='(a,7f11.5,200f11.5)'
      version='DIURNRET   Version 4.14.10  22-Nov-2011   GCT'
      open(lun_rpt,file='diurnret.rpt',status='unknown')
      write(6,*)version
      write(lun_rpt,*)version
c
      zminwas=0.0  ! avoid compiler warnings
      dzmax=0.0  ! avoid compiler warnings
      write(6,117)
 117  format(' Enter Slant Column File (e.g. fai97rat.lav): ',$)
      read(5,'(a)')lavfile
      write(lun_rpt,*) lavfile
      lo=lnbc(lavfile)
c==================================================================
      ieof=0
c  Read slant paths, slant columns, and slant column errors for each spectrum
c      write(*,*)lo,lavfile(:lo-4)//'.mav'
      open(lun_mav,file=lavfile(:lo-4)//'.mav',status='old')
      read(lun_mav,*)   ! read version string
      open(lun_ray,file=lavfile(:lo-4)//'.ray',status='old')
      read(lun_ray,*)nn,nlev_ray
      nlev_ray=nlev_ray-7-ncell  ! First 7 cols are Spectrum Zobs Pres SZA Bend FOV Zmin
      if(nlev_ray.gt.mlay) stop 'Increase parameter MLAY'
      call skiprec(lun_ray,nn-1)
      read(lun_ray,'(a)')ray_string
      open(lun_lav,file=lavfile,status='old')
      read(lun_lav,'(i2,i4,i7,i4)')ncomm,ncollav,nrow,nauxlav

      if(nauxlav.gt.mauxlav) stop 'nauxlav .gt. mauxlav'
      if(ncomm.gt.mcomm) write(6,*)'Warning: Increase parameter MCOMM'
      do i=2,ncomm-1
         read(lun_lav,'(a)') comment(i)
      end do
      read(lun_lav,'(a)')lavlabel
      call substr(lavlabel,lavstarr,2*mgas+nauxlav,nn)
      if(nn.ne.ncollav) then
         write(*,*)'# columns according to file header: ',ncollav
         write(*,*)'# of column labels actually found: ',nn
      stop 'Inconsistent number of columns/labels'
      endif
      ngas=(ncollav-nauxlav)/2
      if(ngas.gt.mgas) then
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
      open(lun_jac,file=lavfile(:lo)//'.jac',status='unknown')
c
c Write .ret file
      ndrec=0
      open(lun_ret,file=lavfile(:lo)//'.ret',status='unknown')
      write(lun_ret,*)ncomm+1,nauxret+2*ngas  
      write(lun_ret,*)version
      do i=2,ncomm-2
         write(lun_ret,'(a)') comment(i)(:lnbc(comment(i)))
      end do
      output_fmt='(i6,f12.6,7f12.5,3(1pe12.4),MMMe12.4)'
      write(output_fmt(29:31),'(i3.3)') nauxret+2*ngas
      write(lun_ret,'(a)') output_fmt
      gas1=lavstarr(nauxlav+1)
      write(lun_ret,'(a)')
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
           read(lun_lav,input_fmt,end=14)
     &     specname(:ls),year(iobs),doty(iobs),hour(iobs),
     &     run(iobs),oblat(iobs),oblon(iobs),zobs(iobs),zmin(iobs),
     &     asza(iobs),azim(iobs),dum,dum,dum,dum,tins(iobs),pins(iobs),
     &     tobs(iobs),pobs(iobs),hobs(iobs),dum,dum,dum,dum,
     &     (slcol(iobs,kgas),slerr(iobs,kgas),kgas=1,ngas)
c        else
c           read(lun_lav,input_fmt,end=14)
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
c        write(*,*) ray_string
        read(ray_string,ray_format) runlab,dum,dum,
     &  obs_asza,dum,dum,dum,(dcell,j=1,ncell),
     &  (sp(iobs,jlay),jlay=1,nlev_ray)
        read(lun_ray,'(a)',end=16) ray_string
        lr=lnbc(runlab)
c        write(*,'(a)')col1//runlab(:lr+3)//runwas(:lr+3)//next_spectrum
c     &//ray_string(:lr)
        if (ray_string(lr-3:lr).eq.'.adf') then
c          write(*,*)'ACE spectra'
           ext=7 ! ACE spectra all have .adf ending, so you need to match .XXX.adf
        else
           ext=3
        endif
        do while(ray_string(lr-ext:lr).eq.runlab(lr-ext:lr))
           read(ray_string,ray_format) runlab,dum,dum,
     &     obs_asza,dum,dum,dum,(dcell,j=1,ncell),
     &     (sp(iobs,jlay),jlay=1,nlev_ray)
           read(lun_ray,'(a)',end=16) ray_string
        end do
16      if(col1.eq.'-' .or. col1.eq.':') go to 2  ! skip bad spectra
c
        if(asza(iobs).ne.obs_asza) write(6,'(a,i4,2f10.3)')
     &   runlab(:24)//ray_string(:24)//
     &  ': Warning: zenith angle mismatch ',iobs,asza(iobs),obs_asza
c
c  Compute Latitude and longitude of tangent points
        if(asza(iobs).le.90.) then
           zeff(iobs)=zobs(iobs)+scht*cos(d2r*asza(iobs))
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
        write(*,*) 'Warning: slant paths are all zero for ',runlab
 76     if(ilay.lt.lowlay) lowlay=ilay  ! lowest atmospheric layer used
        if(ilay.gt. hilay) hilay =ilay  ! highest unused atmospheric layer
 
        if(iobs.le.1 ) start_runlab=runlab
        runwas=runlab
        if(col1.eq.'+') go to 15
      end do   !  iobs=1,mobs
      read(lun_lav,*,end=14)
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
c      write(*,*)'Next Spectrum=',next_spectrum,runlab
c==================================================================
c  Read model (altitudes, temperatures, pressures, & densities).
c  The vmrs, which are in the same file, are not used in the retrieval.
c  They are read merely to carry through the unretrieved gases (or
c  unretrieved levels) to the new output files (.VMR & .MAV)
c
c Only do this if the day has changed. If we are retrieving profiles
c from ascent & sunset spectra acquired on the same day, don't read mavfile.
c      write(*,'(a)') runlab(4:8)//'  '//next_spectrum(4:8)
      do while( runlab(4:8) .ne. next_spectrum(4:8) )
         write(6,*)
         read(lun_mav,'(14x,a)',end=99) next_spectrum 
         write(*,*) next_spectrum
         read(lun_mav,*)nn,ncol_mav,nlev_mav
         nlev_mav=nlev_mav-ncell    ! ignore cell levels
         if(nlev_mav.ne.nlev_ray) stop 'nlev_mav .ne. nlev_ray'
         nvmr=ncol_mav-nauxmav   ! First 4 columns are Z T P D
         if(nvmr.gt.mvmr) then
            write(6,*)' mvmr, nvmr = ',mvmr,nvmr
            write(6,*)' ncol_mav, nauxmav = ',ncol_mav,nauxmav
            stop 'Increase parameter MVMR'
         endif
         call skiprec(lun_mav,nn-2)    ! any remaining comments
         read(lun_mav,'(a)') mavlabel  ! column labels
         call lowercase(mavlabel)
         call substr(mavlabel,mavstarr,mvmr+nauxmav,kcol)
         if( kcol .ne. ncol_mav ) then
            write(*,*)'kcol,ncol_mav=',kcol,ncol_mav
            stop 'column # mismatch in mavfile'
         endif
c
c  Skip cell info
         do i=1,ncell
            read(lun_mav,'(2f7.2,2e11.3,200e11.3)')z(1),t(1),
     &      p(1),d(1),(vmr(j,1),j=1,nvmr)
         end do
c  Read atmospheric level info
         do jlay=1,nlev_mav
            read(lun_mav,'(2f7.2,2e11.3,200e11.3)')z(jlay),t(jlay),
     &      p(jlay),d(jlay),(vmr(j,jlay),j=1,nvmr)
         end do
      end do      !   runlab(4:8) .ne. next_spectrum(4:8) 
c      write(6,'(i4,2x,a)') ocltn,runwas(4:9)
      write(lun_rpt,*) runwas(4:9)
      write(6,*)' Number of different retrieved gases = ',ngas
      write(lun_rpt,*)' Number of different retrieved gases = ',ngas
      write(6,*)' Number of observed spectra used     = ',nobs
      write(lun_rpt,*)' Number of observed spectra used     = ',nobs
      write(6,*)' Lowest atmospheric layer used       = ',lowlay
      write(lun_rpt,*)' Lowest atmospheric layer used       = ',lowlay
      write(6,*)' Highest atmospheric tangent layer   = ',hilay+1
      write(lun_rpt,*)' Highest atmospheric tangent layer   = ',hilay+1
      write(6,*)' Number of retrieved layers          = ',nlev_ray
      write(lun_rpt,*)' Number of retrieved layers          = ',nlev_ray
      write(lun_rpt,*)' GAS       Vbar    Vmda    A_Fit(%) O_Fit(%)'
     &//' C_Fit(%) U_Fit(%)  Frac(%)  SMOO'
c==================================================================
c  Start retrievals, gas by gas
      do 300 kgas=1,ngas  ! Main Loop over retrieval gases
        gas=lavstarr(2*kgas+nauxlav-1)
        if(gas(1:4) .eq. 'tco2') gas=gas(2:)
        molid=0
        do jcol=5,nvmr
           ix=lnbc(mavstarr(jcol))
           ig=index(gas,'_')
           if(ig.gt.0) gas=gas(:ig-1)
           if(gas.eq.mavstarr(jcol)(:ix)) go to 122
           if('1'//gas.eq.mavstarr(jcol)(:ix)) go to 122
        end do
122     molid=jcol-4
c        write(*,*)kgas, gas, molid+nauxmav-1
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  Set up the equation  YY = PD.X  where YY=slcol/slerr & PD=dens*slpath/slerr
c  Also compute mean vmr (VBAR) for profile.
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
        dcfile='/ggg/diurnal/'//runwas(4:8)//'/dcf.'//dcflag//'.'//gas
        inquire(file=dcfile,exist=dcfilexst)
        if(dcfilexst) then
           open(lun_dc,file=dcfile,status='old')
           read(lun_dc,*) nn, nsza
           call skiprec(lun_dc,nn-2)
           read(lun_dc,'(a)')dcflabel
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
           do iobs=1,nobs
              do j=1,nl
                 if(z(j).gt.zmin(iobs)) then
                    if(zmin(iobs).eq.zobs(iobs)) then  !  up looking geometry
                       sza=asin((zobs(iobs)+erad)*sin(d2r*asza(iobs))/
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

c In case of CH3CN use vaa=vbar, otherwise rubbish
        if (gas .eq. 'ch3cn') then
           do j=1,nl         
              vaa(j) = vbar   
           end do            
        endif

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
         smfact=0.015  ! MATMOS
c         smfact=0.0625*cos(mtlat)+0.01*(1-cos(mtlat))**2  ! MkIV balloon
c         smfact=0.020  !  for ACE tropical data
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

c  In case of CH3CN retrieval:
c  Augment matrix of partials (PD) with NLAB vbar above the balloon
c  Augment measurement vector (YY) with NLAB a priori constraints above the balloon
        if (gas .eq. 'ch3cn') then
          aprierr=1e-11
          nlab=nlev_ray-hilay-1 ! no. of levels above balloon
          neq=nobs+nl+nlab ! no. of equations
          do ilay=1,nlab
             do jlay=1,nl
                pd1(nobs+nl+ilay,jlay)=0.0
                if(ilay+hilay+1 .eq. jlay+lowlay-1) then
                   pd1(nobs+nl+ilay,jlay)=vbar/aprierr
                endif
             end do
c            write(*,*)'vmr',vmr(molid+nauxmav-1,ilay+hilay+1)
             pd1(nobs+nl+ilay,nl+1)=
     &       vmr(molid+nauxmav-1,ilay+hilay+1)/aprierr
c            write(*,*)'pd1right',pd1(nobs+nl+ilay,nl+1)
          end do  ! ilay=1,nlab
        endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c   Duplicate entire array of partial differentials
        write(lun_jac,'(4i5)')ocltn,kgas,nobs,nl
        do iobs=1,neq
           do ilay=1,nl+1
              pd2(iobs,ilay)=pd1(iobs,ilay)
           end do
           write(lun_jac,'(200(1pe12.4))')(pd1(iobs,jlay),jlay=1,nl+1)
        end do
c
c  Perform actual retrievals; the first constrained, the second unconstrained
        call nnls(pd1,mm,neq,nl,pd1(1,nl+1),vcon(1,kgas),
     &  rnorm,wk,zz,jpvt1,mode)
        call snnls(pd2,mm,neq,nl,vunc(1,kgas),np,jpvt2,wk,0)
        rnorm=slsfin(pd2,mm,neq,nl,np,jpvt2)
        do ilay=1,nl
          xx=(abs(vaa(ilay))+vbar)/2
          vcon(ilay,kgas)=xx*vcon(ilay,kgas)
          vunc(ilay,kgas)=xx*vunc(ilay,kgas)
          ev(ilay,kgas)=
     &    xx*dsqrt(dsdot(nl,pd2(ilay,1),mm,pd2(ilay,1),mm))
        end do

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  Recalculate slant columns using original and both constrained (VCON) &
c  unconstrained (VUN)  retrieved vmr profiles.
c  Compute the "goodness of fit" to the measured slant columns.
c  Compute the factor (FRAC) by which the deviations of the measurements
c  from the re-calculated slant columns exceeds their original error estimates.
        tverr=0.0
        ta=0.0
        to=0.0
        tc=0.0
        tu=0.0
        tt=0.0
        do iobs=1,nobs
          vaacol=0.0
          orgcol=0.0
          concol=0.0
          unccol=0.0
          do ilay=1,nl
            xx=d(ilay+lowlay-1)*km2cm*sp(iobs,ilay+lowlay-1)
            vaacol=vaacol+xx*vaa(ilay)
            unccol=unccol+xx*vunc(ilay,kgas)
            concol=concol+xx*vcon(ilay,kgas)
            orgcol=orgcol+xx*vmr(molid,ilay+lowlay-1)
          end do
c          write(6,*)orgcol,slcol(iobs,kgas),concol,unccol
          tverr=tverr+(vaacol/slerr(iobs,kgas))**2
          ta=ta+((slcol(iobs,kgas)-vaacol)/slerr(iobs,kgas))**2
          to=to+((slcol(iobs,kgas)-orgcol)/slerr(iobs,kgas))**2
          tc=tc+((slcol(iobs,kgas)-concol)/slerr(iobs,kgas))**2
          tu=tu+((slcol(iobs,kgas)-unccol)/slerr(iobs,kgas))**2
          tt=tt+ (slcol(iobs,kgas)/slerr(iobs,kgas))**2
c          write(6,*)ta,to,tc,tu,tt
        end do
        afit=sqrt(ta/tt)
        ofit=sqrt(to/tt)
        cfit=sqrt(tc/tt)
        ufit=sqrt(tu/tt)
        frac=sqrt(tu)/nobs
        verr=vbar/sqrt(tverr)
      write(lun_rpt,'(a10,1p2e9.2,0p5f8.2,0pf7.1,a2)')
     &  lavstarr(2*kgas+nauxlav-1),vbar,
     &  verr,100*afit,100*ofit,100.*cfit,100*ufit,100*frac,smoo,dcflag
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  Place retrieved profiles in vmr array, except for CO2 and N2. Don't change them.
        if(gas(:lnbc(gas)).ne.'co2' .and. gas(:lnbc(gas)).ne.'n2') then
          vmin=0.001*vbar
          do i=1,nl
            if(vcon(i,kgas).le.vmin) then
              vmr(molid,i+lowlay-1)=vmin
            else
              vmr(molid,i+lowlay-1)=vcon(i,kgas)
            endif
          end do
        endif
 300  continue    ! kgas=1,ngas   Main Loop over retrieval gases
c==================================================================
c  Create matrix of retrieved vmrs in a format which can be plotted by GPLOTEM
c  Plot them in chronological order.
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
          aszai=aszai+wt*asza(iobs)
          aszaz=aszaz+wt*asza(iobs)*dz
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
        write(lun_ret,output_fmt)
     &  ocltn,yeari,dotyi,houri,tlati,tloni,aszai,
     &  z(ilay+lowlay-1),t(ilay+lowlay-1),1013.25*p(ilay+lowlay-1),
     &  d(ilay+lowlay-1),
     &  t(ilay+lowlay-1)*p(ilay+lowlay-1)**(-.286),
     &  (vunc(ilay,kgas),ev(ilay,kgas),kgas=1,ngas)
c
c  PD1 contains the constrained retrieved vmrs 
c  PD2 contains the unconstrained retrieved vmrs
c  SP  contains their uncertainties.
      end do   ! ilay=2,hilay-lowlay+1+3
c=======================================================================
      if(ieof.gt.0) go to 99
      end do  ! ocltn=1,999999
99    close(lun_ray)
      close(lun_mav)
      close(lun_lav)
      close(lun_rpt)
      close(lun_jac)
      close(lun_ret)
      end
