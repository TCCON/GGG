c  Program to convert old-format runlogs into new runlogs.
c
      implicit none

      real*8 
     & oblat,            ! observation latitude (deg)
     & oblon,            ! observation longitude (deg)
     & obalt,            ! observation altitude (km)
     & asza,             ! astronomical solar zenith angle (unrefracted)
     & azim,             ! solar azimuth angle (unrefracted)
     & osds,             ! object-source doppler stretch
     & opd,              ! Optical path difference (cm) of interferogram
     & graw,             ! spacing of raw spectrum (cm-1) from GETINFO
     & zpdtim,           ! Time of ZPD (UT hours)
     & zenoff,           ! Zenith angle pointing offset (deg)
     & fovi,             ! Internal angular diameter of FOV (radians)
     & fovo,             ! External angular diameter of FOV (radians)
     & amal,             ! angular misalignment of interferometer (radians)
     & zoff,             ! Zero level offset (dimensionless fraction)
     & snr,              ! Signal-To-Noise Ratio (dimensionless)
     & tins,             ! Inside temperature (C)
     & pins,             ! Inside pressure (mbar)
     & hins,             ! Inside humidity (%)
     & tout,             ! Outside temperature (C)
     & pout,             ! Outside pressure (mbar)
     & hout,             ! Outside humidity (%)
     & wspd,             ! Wind Speed (m/s)
     & wdir,             ! Wind Direction (deg.)
     & sia,              ! Solar Intensity (Average)
     & fvsi,             ! Fractional Variation in Solar Intensity
     & aipl,             ! Airmass-Independent Path Length (km)
     & lasf,             ! Laser Frequency (e.g. 15798 cm-1)
     & wavtkr            ! Suntracker frequency (active tracking)

      character
     & col1*1,           ! first column
     & specname*(57), ! spectrum name
     & apf*2,            ! apodization function (e.g. BX N2, etc)
     & version*50,
     & data_fmt_read_rl*256,
     & data_fmt_write_rl*256,
     & col_headers_rl*320,
     & infile*24,
     & outfile*24

      integer*4 ispec,
     & lunw_rlg,         ! Logical Unit Number
     & lunr_rlg,         ! Logical Unit Number
     & istat,            ! status flag (0=success, 1=EOF)
     & iyr,              ! year
     & iset,             ! day of year
     & ifirst,           ! index of first spectral point in disk file
     & ilast,            ! index of last spectral point in disk file
     & possp,            ! Length of attached header in bytes
     & bytepw            ! Bytes per data word, usually 2 (I*2) or 4 (R*4)

      parameter (lunr_rlg=14,lunw_rlg=16)

      version =' convert_runlog   Version 0.11   2016-04-02   GCT '
      col1=' '
      outfile='convert_runlog.out'
      if (iargc() == 0) then
         write(*,'(a)')'Name of old runlog:'
         read(*,'(a)') infile
      elseif (iargc() == 1) then
         call getarg(1, infile)
      else
         stop 'Usage: $gggpath/bin/convert_runlog old_runlog'
      endif
      open(lunr_rlg,file=infile,status='old')
      call read_runlog_header(lunr_rlg, data_fmt_read_rl,col_headers_rl)
c      read(lunr_rlg,*)nlhead,ncol
c      do j=2,nlhead
c         read(lunr_rlg,*)
c      enddo

      open(lunw_rlg,file=outfile,status='unknown')
      call write_runlog_header(lunw_rlg,version,data_fmt_write_rl)

c      write(lunw_rlg,*) ' 3   36'
c      write(lunw_rlg,*) 'convert_runlog'
c      write(lunw_rlg,'(a)')' Spectrum_File_Name                     Year'//
c     &  '  Day   Hour   oblat    oblon   obalt    ASZA    POFF   AZIM'//
c     &  '  OSDS    OPD   FOVI  FOVO  AMAL   IFIRST    ILAST'//
c     &  '    DELTA_NU    POINTER  BPW ZOFF SNR  APF tins  pins'//
c     &  '  hins   tout   pout  hout  sia    fvsi   wspd  wdir'//
c     &  '  lasf    wavtkr  aipl'

      do ispec=1,9999999

          call read_runlog_data_record(lunr_rlg,data_fmt_read_rl,
     &    col1,specname,iyr,iset,zpdtim,
     &    oblat,oblon,obalt,asza,zenoff,azim,osds,opd,fovi,fovo,amal,
     &    ifirst,ilast,graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     &    tout,pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)

          if (istat.ne.0) exit

          call write_runlog_data_record(lunw_rlg,data_fmt_write_rl,
     &    col1,specname,iyr,iset,zpdtim,
     &    oblat,oblon,obalt,asza,zenoff,azim,osds,opd,fovi,fovo,amal,
     &    ifirst,ilast,graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     &    tout,pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)
      end do  !  ispec=1,9999999

      close(lunr_rlg)
      close(lunw_rlg)
      stop
      end
