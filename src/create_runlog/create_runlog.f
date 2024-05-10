c  Program to create a GGG-compatible runlog from the sunrun file.
c  which is created by running create_sunrun_from_xxxx.
c
      implicit none
      include "../gfit/ggg_int_params.f"

      integer*4
     & bytepw,   !Number of bytes per data word (=2 for AT,M4; =4 for OP,GR)
     & ifirst,   !The index of the first spectral point on disk
     & ilast,    !The index of the last spectral point on disk
     & possp,    !Number of bytes before IFIRST'th point (i.e. header length)
     & object    !Heavenly object (Moon=1, Sun=2)

      integer*4
     & iy,im,id,jd,
     & jj,
     & i,idum,ios,
     & lnbc,lloc,ispe,ifmin,ifmax,
     & nlhead,ncol,
     & lr,lrt,
     & doy,
     & one,
     & lunr,lunw_rlg,lst,lunw_lse,lunw_nts
      parameter (lunr=15,lunw_rlg=16,lunw_lse=17,one=1,lunw_nts=62)
c
      real*8 amal,fovi,fovo,gmt,tins,pins,hins,tout,pout,hout,
     & obalt,asza,azim,snr,wavtkr,zoff,zpoff,oblat,oblon,opd,
     & lasf,fmin,fmax,fsf,delwav,tcorr,sia,fvsi,wspd,wdir,aipl,
     & tel_mag,eorv,ervc,osds,tplat,tplon,tpalt,site_solar_noon,
     & lse,lsu,lsf,dip,mvd,r8was,r8year,r8ydiff
c
      character spfmt*2, logfile*40,
     & outfile*(mfilepath),path*(mpath),dplist*(mfilepath),
     & lsefile*(mfilepath)

      character
     & col1*1,                    !first column of runlog record
     & apf*2,                     !apodization function (e.g. BX N2, etc)
     & dl*1,                      !forward or backward slash
     & ext*3,                     !geometry ['air','bal','gnd','lab',orb','syn']
     & gggdir*(mpath),            !ggg directory path (GGGPATH?)
     & specname*(nchar),          !spectrum name
     & data_fmt_wrlg*256,
     & lsefile_format*100,
     & version*64                 !current program version


      idum=mauxcol ! Avoid compiler warning (unused parameter)
      idum=mcolvav ! Avoid compiler warning (unused parameter)
      idum=mgas    ! Avoid compiler warning (unused parameter)
      idum=mlev    ! Avoid compiler warning (unused parameter)
      idum=mrow_qc ! Avoid compiler warning (unused parameter)
      idum=mspeci  ! Avoid compiler warning (unused parameter)
      idum=mvmode  ! Avoid compiler warning (unused parameter)
      idum=ncell   ! Avoid compiler warning (unused parameter)

c
      version=' CREATE_RUNLOG    Version 8.78     2019-08-22    GCT'
      write(*,'(a)') version
      col1=' '
      iy=0
      im=0
      id=30
      amal=0.0d0
      r8was=-9999999.9d0


      call get_ggg_environment(gggdir, dl)
      lrt=lnbc(gggdir)       !Length of gggdir

      lr=0
      do while(lr.eq.0)
         if (iargc() == 0) then
            write(6,'(a)') 'Enter sunrun (e.g. fts89avg.gop): '
            read(*,'(a)') logfile
         elseif (iargc() == 1) then
            call getarg(1, logfile)
         else
            stop 'Usage: $gggpath/bin/create_runlog sunrun'
         endif
         lr=lnbc(logfile)
      end do
      ext(1:1)=logfile(lr-2:lr-2)
      if(ext(1:1).eq.'a') ext(1:3)='air'
      if(ext(1:1).eq.'b') ext(1:3)='bal'
      if(ext(1:1).eq.'g') ext(1:3)='gnd'
      if(ext(1:1).eq.'l') ext(1:3)='lab'
      if(ext(1:1).eq.'o') ext(1:3)='orb'
      if(ext(1:1).eq.'s') ext(1:3)='syn'

c
      spfmt=logfile(lr-1:lr)
      open(lunr,file=gggdir(:lrt)//'sunruns'//dl//ext//dl//logfile,
     $status='old')
      read(lunr,*) nlhead, ncol
      do i=2,nlhead
         read(lunr,*)
      end do
      outfile=gggdir(:lrt)//'runlogs'//dl//ext//dl//logfile(:lr-2)//'rl'
      lsefile=gggdir(:lrt)//'lse'//dl//ext//dl//logfile(:lr-3)//'lse'
c     write(*,*)'lsefile=',lsefile
c
      open(lunw_rlg,file=outfile,status='unknown')
      open(lunw_lse,file=lsefile,status='unknown')
      open(lunw_nts,file='create_runlog.nts',status='unknown')
      lsefile_format = 
     & '(1x,a99,i5,i4,f9.4,i3,1pe12.4,1pe12.4,1pe12.4,'//
     & '1pe12.4,0p,1x,f9.4)'
      write(lsefile_format(6:7),'(i2.2)') nchar
      write(lunw_lse,'(a)')' 3  10'
      write(lunw_lse,'(a)')' '//lsefile_format
      write(lunw_lse,'(a)')
     &' Specname                                                 '//
     &' year doy  hour    LST LSE         LSU         LSF        '//
     &' DIP          MVD'
      call write_runlog_header(lunw_rlg,version,data_fmt_wrlg)
c
      do ispe=1,9999999  !---------Main loop over spectra----------
         call read_sunrun(lunr,col1,specname,object,tcorr,oblat,
     &   oblon,obalt,tins,pins,hins,tout,pout,hout,sia,fvsi,
     &   wspd,wdir,fmin,fmax,fsf,lasf,wavtkr,aipl,tel_mag,ios)
c        write(*,*)ispe,specname,ios
         if(ios.ne.0) go to 99
c
c  find the spectral file, return the PATH to the spectrum
         dplist=gggdir(:lrt)//'config'//dl//'data_part.lst'
         call gindfile(dplist,specname,path)
         if(lnbc(path).le.0) then
            write(6,*) ' Not Found: "'//specname//'"'
            stop
         endif

c         write(*,*)'create_runlog: ',path,specname,iy,im,id,gmt
         call rdsphead(spfmt,specname,path(:lloc(path,dl)),
     &   ifirst,ilast,possp,
     &   bytepw,apf,delwav,opd,fovi,snr,
     &   iy,im,id,gmt,lasf,wavtkr,lst,lse,lsu,lsf,dip,mvd)

         fovo=fovi/tel_mag
c         write(*,*) path,ifirst,ilast,possp,bytepw,apf,delwav

c  Fudge the starting index to prevent use of noise regions.
c  Multiply by delwav to allow for reverse-ordered spectra
c  in which IFIRST, ILAST & DELWAV are all negative.
         ifmin=nint(fmin/delwav)
         if(ifmin*delwav.gt.ifirst*delwav) then
            possp=possp+iabs(bytepw)*(ifmin-ifirst)
            ifirst=ifmin
         endif
c
c  Fudge the ending index to prevent use of noise regions.
         ifmax=nint(fmax/delwav)
c         write(*,*)ifmax,ilast,fmax,ilast*delwav
         if(ifmax*delwav.lt.ilast*delwav) then
            ilast=ifmax
         endif
c
         delwav=delwav*fsf
         gmt=gmt+tcorr/3600.0d0
c
c  Calculate Day of Year (DOY)
         call julian(iy,im,id,jd)
         call julian(iy,one,one,jj)
         doy=jd-jj+1
c
c  We want all observations to have the day number of local noon,
c  even though the UT day number may change. The following code
c  makes all spectra acquired on the same local day have the same
c  day number. To achieve this, the UT time must be allowed to go
c  -ve or to exceed 24.
         site_solar_noon=12-oblon/15  !  ut time of local noon
         if((gmt-site_solar_noon).gt.+12) then
            gmt=gmt-24.0
            doy=doy+1
         elseif((gmt-site_solar_noon).lt.-12) then
            gmt=gmt+24.0
            doy=doy-1
         endif
c         write(*,*)path(:lnbc(path))
c
         if(object.gt.0) then
            call zenaz(object,oblat,oblon,obalt,iy,one,doy,
     &      gmt/24.d0,asza,azim,eorv,ervc,tplat,tplon,tpalt)
            osds=1.e+06*(eorv+ervc)/3.d+8  ! Observer-Sun Doppler Stretch (ppm)
c            write(*,*)'osds=',osds,eorv,ervc
         else
            asza=0.0
            azim=0.0
            osds=0.0
         endif

         zoff=0.0
         zpoff=0.0
         call write_runlog_data_record(lunw_rlg,data_fmt_wrlg,
     &   col1,specname,iy,doy,gmt,oblat,oblon,obalt,asza,zpoff,azim,
     &   osds,opd,fovi,fovo,amal,ifirst,ilast,
     &   delwav,possp,bytepw,zoff,snr,apf,tins,pins,hins,tout,
     &   pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,ios)

c  Create report to identify non-chronological times in the runlog
         r8year=iy+(doy+gmt/24.0d0)/366.0d0
         r8ydiff=r8year-r8was
c  Report negative time-steps in the runlogs times of more than 1.5min.
         if( r8ydiff .lt. -0.000003d0) then
           write(lunw_nts,'(a,a,2f12.6)')
     &  '  Negative time step (runlog not chronologically sorted): ',
     &     specname,r8was,r8year
         endif
         r8was=r8year


         write(lunw_lse,lsefile_format)
     &   specname,iy,doy,gmt,lst,lse,lsu,lsf,dip,mvd
         if(mod(ispe,1000).eq.0) write(*,*) ispe
c
c         write(*,*)'---------------------------------------------------'
      end do ! -------------Main loop over spectra----------------
c
 99   close(lunr)
      close(lunw_rlg)
      close(lunw_lse)
      close(lunw_nts)

      write(*,*) ispe-1, ' spectra found'
      stop
      end
