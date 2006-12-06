      subroutine write_runlog(lun,col1,runlab,iyr,iset,zpdtim,
     & oblat,oblon,obalt,asza,zenoff,opd,fovi,fovo,amal,ifirst,
     & ilast,graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     & tout,pout,hout,lasf,wavtkr,sia,sis,aipl,istat)
c
c  Writes a single record to a runlog file.
c  The file must already have been opened.
c  The goal of this subroutine is to hide all the of the code
c  that depends on the runlog format into pair of subroutines
c  (writrunlog & readrunlog).
c
c  This has two advantages:
c  1) It simplifies the calling programs
c  2) It means that if the runlog format is changed
c  in the future, only write_runlog &read_runlog needs to be changed,
c  not the dozen main programs that read the runlog.
c
c  Input:
c    lun   Logical Unit number of file to be read
c    everything else

c  Outputs:
c    istat

      integer*4
     & lun,              ! Logical unit number
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
     & obalt,             ! observation altitude (km)
     & asza,             ! astronomical solar zenith angle (unrefracted)
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
     & sia,              ! Solar Intensity (Average)
     & sis,              ! Solar Intensity (SD)
     & aipl,             ! Airmass-Independent Path Length (km)
     & lasf,             ! Laser Frequency (e.g. 15798 cm-1)
     & wavtkr            ! Suntracker frequency (active tracking)

      character
     & col1*1,           ! first column
     & runlab*21,        ! spectrum name
     & apf*2             ! apodization function (e.g. BX N2, etc)

      write(lun,34,err=99) col1,runlab,iyr,iset,zpdtim,oblat,oblon,
     &obalt,asza,zenoff,opd,fovi,fovo,amal,ifirst,ilast,graw,possp,
     &bytepw,zoff,nint(snr),apf,tins,pins,hins,
     &tout,pout,hout,lasf,wavtkr,sia,sis,aipl
 34   format(a1,a21,1x,2i4,f8.4,f8.3,f9.3,2f8.3,1x,f6.4,f7.2,3(1x,f5.4),
     &2i8,1x,f14.11,i8,i3,1x,f5.3,i5,1x,a2,2(f6.1,f8.2,f5.1),f10.3,f7.0,
     & 2f6.1,f7.3)
      istat=0
      return
99    istat=1
      return
      end
