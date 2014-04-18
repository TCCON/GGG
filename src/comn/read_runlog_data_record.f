      subroutine read_runlog_data_record(lunr_rlg,data_fmt_read_rl,
     & col1,specname,iyr,iset,zpdtim,
     & oblat,oblon,obalt,asza,zenoff,azim,osds,opd,fovi,fovo,amal,
     & ifirst,ilast,graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     & tout,pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)
c
c  Reads a single record from the runlog file.
c  File must already have been opened.
c  The purpose of this subroutine is to hide all the
c  of the code that depends on the runlog format into
c  a single subroutine. This has two advantages:
c  1) It simplifies the calling programs
c  2) It means that if the runlog Format is changed
c  in the future, only read_runlog subroutine needs to be changed,
c  not the dozen main programs that read the runlog.
c
c
c  Input:
c    lunr_rlg          I*4   Logical Unit number of file to be read
c    data_fmt_read_rl  C*(*) Runlog data format

c  Outputs:
c    everything else

      implicit none
      integer*4 lnbc,isnr,
     & lunr_rlg,         ! Logical unit number
     & istat,            ! status flag (0=success, 1=EOF)
     & iyr,              ! year 
     & iset,             ! day of year
     & ifirst,           ! index of first spectral point in disk file
     & ilast,            ! index of last spectral point in disk file
     & possp,            ! Length of attached header in bytes
     & bytepw            ! Bytes per data word, usually 2 (I*2) or 4 (R*4)

      real*8
     & oblat,            ! observation latitude (deg).
     & oblon,            ! observation longitude (deg).
     & obalt,            ! observation altitude (km)
     & asza,             ! astronomical solar zenith angle (unrefracted)
     & opd,              ! Optical path difference (cm) of interferogram
     & graw,             ! spacing of raw spectrum (cm-1) from GETINFO
     & zpdtim,           ! Time of ZPD (UT hours)
     & zenoff,           ! Zenith angle pointing offset (deg)
     & azim,             ! Solar azimuth angle
     & osds,             ! Observer-Sun Doppler Shift (ppm)
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
     & sia,              ! Solar itensity
     & fvsi,             ! Fractional Variation in Solar Intensity
     & wspd,             ! Wind Speed (m/s)
     & wdir,             ! Wind Direction (deg)
     & aipl,             ! Airmass-Independent Path Length (km)
     & lasf,             ! Laser Frequency (e.g. 15798 cm-1)
c     & zerr,             ! Extra parameter in ATMOS tab-delimited runlogs
c     & scalf,            ! Extra parameter in ATMOS tab-delimited runlogs
     & wavtkr            ! suntracker frequency (active tracking)

      character
     & col1*1,           ! first column of runlog record
c     & record*400,       ! runlog record
     & data_fmt_read_rl*256,
     & specname*(*),     ! spectrum name
     & apf*2             ! apodization function (e.g. BX N2, etc)

c       read(lunr_rlg,'(a1,a)',end=99) col1,record
c1      read(lunr_rlg,'(a1,a)',end=99) col1,record
c      if( col1.eq.':') go to 1
c      if(col1.ne.'-' .and. col1.ne.'+' .and. col1.ne.' ') then
c         record=col1//record   ! Runlog is the old format
c      endif
c      lr=lnbc(record)
c      write(*,*) record
c      write(*,*)'read_runlog: lr= ', lr

c   Format was successfully read from runlog header, so use it.
      if(lnbc(data_fmt_read_rl).gt.0) then  
        read(lunr_rlg,data_fmt_read_rl,end=98,err=99)col1,specname,
     &  iyr,iset,zpdtim,oblat,oblon,
     &  obalt,asza,zenoff,azim,osds,opd,fovi,fovo,amal,ifirst,ilast,
     &  graw,possp,bytepw,zoff,isnr,apf,tins,pins,hins,tout,pout,hout,
     &  sia,fvsi,wspd,wdir,lasf,wavtkr,aipl
        snr=float(isnr)
        istat=0
      else                 !
        istat=1
        write(*,*) 'read_runlog_data_record: Missing runlog format'
      endif
      return
98    istat=2
      return
99    istat=3
      return
      end
