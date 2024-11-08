      subroutine write_runlog_data_record(lunw_rlg,data_fmt_wrlg,
     & col1,specname,iyr,iset,zpdtim,
     & oblat,oblon,obalt,asza,zenoff,azim,osds,opd,fovi,fovo,amal,
     & ifirst,ilast,graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     & tout,pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)
c
c  Writes one data record to the runlog file (lunw_rlg),
c  which must already have been opened and the header lines
c  written (by calling write_runlog_header).
c  The goal of this subroutine is to hide all the of the code
c  that depends on the runlog format into pair of subroutines
c  (write_runlog & read_runlog).
c
c  This has two advantages:
c  1) It simplifies the calling programs
c  2) It means that if the runlog format is changed in the future,
c     only write_runlog & read_runlog subroutines need be changed,
c     not the dozen programs that use the runlog.
c
c  Input:
c    lunw_rlg       Logical Unit number of runlog file
c    data_fmt_wrlg  Format statement (defined in write_runlog_header)
c    everything else

c  Outputs:
c    istat

c  Typical Usage:
c    open(lunw_rlg,file='runlog_name')
c    call write_runlog_header(lunw_rlg,header,data_fmt_wrlg)
c    do ispec=1,nspec
c       call write_runlog_data_record(lunw_rlg,data_fmt_wrlg, ....)
c    end do
c    close(lunw_rlg)

      implicit none

      integer*4
     & ksnr,
     & lunw_rlg,              ! Logical unit number
     & istat,            ! status flag (0=success, 1=EOF)
     & iyr,              ! year 
     & iset,             ! day of year
     & ifirst,           ! index of first spectral point in disk file
     & ilast,            ! index of last spectral point in disk file
     & possp,            ! Length of attached header in bytes
     & bytepw            ! Bytes per data word, usually 2 (I*2) or 4 (R*4)

      real*8  wlimit,
     & oblat,            ! observation latitude (deg).
     & oblon,            ! observation longitude (deg).
     & obalt,            ! observation altitude (km)
     & asza,             ! astronomical solar zenith angle (unrefracted)
     & azim,             ! solar azimuth angle
     & osds,             ! Observer-Sun Doppler stretch (ppm)
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
     & sia,              ! Solar Intensity Average (arbitrary units)
     & fvsi,             ! Fractional Variation in Solar Intensity
     & wspd,             ! Wind Speed (m/s)
     & wdir,             ! Wind Direction (deg)
     & aipl,             ! Airmass-Independent Path Length (km)
     & lasf,             ! Laser Frequency (e.g. 15798 cm-1)
     & wavtkr            ! Suntracker frequency (active tracking)

      character
     & col1*1,           ! first column
     & specname*(*),     ! spectrum name
     & data_fmt_wrlg*(*), ! runlog data format
     & apf*2             ! apodization function (e.g. BX N2, etc)

      ksnr=nint(snr)
      if(ksnr.lt.-999) ksnr=-999
      if(ksnr.gt.9999) ksnr=9999
      write(lunw_rlg,data_fmt_wrlg,err=99) col1,specname,
     & iyr,iset,zpdtim,
     & oblat,oblon,obalt,asza,zenoff,azim,osds,wlimit(opd,'f7.2'),
     & fovi,fovo,amal,ifirst,ilast,graw,possp,
     & bytepw,zoff,ksnr,apf,tins,pins,hins,tout,pout,hout,
     & sia,fvsi,wspd,wdir,lasf,wavtkr,aipl
      istat=0
      return
99    istat=1
      return
      end
