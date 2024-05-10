c  Program for writing zenith pointing offsets, previously calculated by
c  ZENANG, to the runlog files (both HgCd and InSb).
c
      implicit none
      include "../gfit/ggg_int_params.f"

      integer*4 lunr_rlg,lunw_rlg,lunz,lnbc,lr,istat,ndet,jdet,iflag,
     & ifirst,ilast,iyr,iset,possp,bytepw,ispe,ntot,idum
      parameter (lunr_rlg=14,lunw_rlg=15,lunz=16)
c
      real*8 year,yearz,aszaz,poutz,totcon,del,
     & td,td2,
     & totwas,totzen,
     & totwas2,totzen2

      real*8
     & oblat,           ! observation latitude (deg).
     & oblon,           ! observation longitude (deg).
     & obalt,            ! observation altitude (km)
     & zpdtim,           ! Time of ZPD (UT hours)
     & asza,             ! astronomical solar zenith angle (unrefracted)
     & zenoff,           ! zenith pointing offset
     & graw,             ! spacing of raw spectrum (cm-1) 
     & opd,              ! Optical path difference (cm) of interferogram
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
c    & sis,              ! Solar Intensity (SD)
     & fvsi,             ! Fractional Variation in Solar Intensity
     & wspd,             ! Wind Speed (m/s)
     & wdir,             ! Wind Direction (deg)
     & azim,             ! Solar Azimuth Angle
     & osds,             ! Observer-Sun Doppler Stretch (ppm)
     & aipl,             ! Airmass-Independent Path Length (km)
     & lasf,             ! laser frequency (e.g. 15798.03 cm-1)
     & wavtkr            ! suntracker operating frequency (e.g. 9900 cm-1)

c
      character
c     & title*240,
     & data_fmt_read_rl*256,
     & data_fmt_write_rl*256,
     & version*52,
     & col_labels_rl*320,
     & col1*1,                    !first column of runlog record
     & apf*2,                     !apodization function (e.g. BX N2, etc)
     & dl*1,
     & gggdir*(mpath),            !ggg directory path (GGGPATH?)
     & specname*(nchar),          !spectrum name
     & rlgfile*120                !name of runlog file

      idum=mfilepath ! Avoid compiler warning (unused parameter)
      idum=mauxcol  ! Avoid compiler warning (unused parameter)
      idum=mcolvav  ! Avoid compiler warning (unused parameter)
      idum=mgas     ! Avoid compiler warning (unused parameter)
      idum=mlev     ! Avoid compiler warning (unused parameter)
      idum=mrow_qc  ! Avoid compiler warning (unused parameter)
      idum=mspeci   ! Avoid compiler warning (unused parameter)
      idum=mvmode   ! Avoid compiler warning (unused parameter)
      idum=ncell    ! Avoid compiler warning (unused parameter)
      idum=nchar    ! Avoid compiler warning (unused parameter)
c
      iflag=0  ! avoid compiler warning (may be used uninitialized)
      write(6,*) version
      version=' ZENCORR     Version 1.53     2016-04-02      GCT '
      call get_ggg_environment(gggdir, dl)

      if (iargc() == 0) then
         write(6,'(a)') 'Enter runlog (e.g. flt93avg.brl) '
         read(5,'(a)') rlgfile
      elseif (iargc() == 1) then
         call getarg(1, rlgfile)
      else
         stop 'Usage: $gggpath/bin/zencorr runlogname'
      endif
      lr=lnbc(rlgfile)
      if(lr.le.0) stop ' Empty runlog name'
c
c  Open the original runlog
      if(rlgfile(lr-2:lr-2).eq.'b') then
         open(lunr_rlg,file=gggdir(:lnbc(gggdir))//'runlogs'//dl//'bal'
     &    //dl//rlgfile, status='old')
         ndet=2
      elseif (rlgfile(lr-2:lr-2).eq.'o') then
         open(lunr_rlg,file=gggdir(:lnbc(gggdir))//'runlogs'//dl//'orb'
     &    //dl//rlgfile, status='old')
         ndet=1
      else
         stop 'unknown extention'
      endif

      call read_runlog_header(lunr_rlg,data_fmt_read_rl,col_labels_rl)
c      read(lunr_rlg,*) nhl,ncol
c
c  Open the new runlog
      open(lunw_rlg,file='new_rl.out',status='unknown')
      call write_runlog_header(lunw_rlg,version,data_fmt_write_rl)
c      write(lunw_rlg,*) nhl,ncol
c      
c      do i=2,nhl
c         read(lunr_rlg,'(a)') title
c         write(lunw_rlg,'(a)') title(:lnbc(title))
c      end do
c
c  Also open the ZENANG.OUT file containing the pointing offsets
      open(lunz,file='zenang.out',status='old')

      ntot=0
      td=0.0d0
      td2=0.0d0
      totwas=0.0d0
      totwas2=0.0d0
      totzen=0.0d0
      totzen2=0.0d0
      read(lunz,*) yearz,totcon,poutz,aszaz,del
      do ispe=1,999999
         do jdet=1,ndet
         call read_runlog_data_record(lunr_rlg,data_fmt_read_rl,
     &   col1,specname,iyr,iset,zpdtim,
     &   oblat,oblon,obalt,asza,zenoff,azim,osds,
     &   opd,fovi,fovo,amal,ifirst,ilast,
     &   graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     &   tout,pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)
         if(istat.gt.0) go to 99
         if(zenoff.lt.-.9999)zenoff=-.9999
         if(zenoff.gt.9.9999)zenoff=9.9999
c        write(*,*)'specname ',specname,totcon,aszaz,asza+zenoff
c        write(*,*)specname,istat,totcon,aszaz,asza+zenoff
         year=iyr+(iset+zpdtim/24)/366.00
         write(*,*) ispe,jdet,year,yearz,asza+zenoff,aszaz
         if(abs(abs(year)-abs(yearz)).lt.0.000001 .and.
     &   abs(asza+zenoff-aszaz).lt.0.0001) then
            iflag=0
            ntot=ntot+1
            td=td+del
            td2=td2+del**2
            totwas=totwas+zenoff
            totwas2=totwas2+zenoff**2
            if(col1.ne.'-') zenoff=zenoff+del
            totzen=totzen+zenoff
            totzen2=totzen2+zenoff**2
         else
            iflag=1
            write(*,*) 'Spectrum mismatch',abs(year),abs(yearz)
            write(*,*) 'Spectrum mismatch',asza+zenoff,aszaz
         endif
         call write_runlog_data_record(lunw_rlg,data_fmt_write_rl,
     &   col1,specname,iyr,iset,zpdtim,
     &   oblat,oblon,obalt,asza,zenoff,azim,osds,opd,fovi,fovo,amal,
     &   ifirst,ilast,graw,possp,bytepw,zoff,snr,apf,tins,pins,
     &   hins,tout,pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)
         end do  ! jdet=1,ndet
         if(iflag.eq.0) read(lunz,*,end=67) yearz,totcon,poutz,aszaz,del
67       continue
      end do  ! ispe=1,99999
c
 99   close(lunz)
      close(lunr_rlg)
      close(lunw_rlg)
      write(*,*)'NTOT = ',ntot,ispe-1
      write(*,*)'Average ZENOFF was = ',totwas/ntot,' +/- ',
     & sqrt(ntot*totwas2-totwas**2)/ntot
      write(*,*)'Average ZENOFF now = ',totzen/ntot,' +/- ',
     & sqrt(ntot*totzen2-totzen**2)/ntot
      write(*,*)'Average DEL now = ',td/ntot,' +/- ',
     & sqrt(ntot*td2-td**2)/ntot
      stop
      end
