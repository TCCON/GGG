c  Program to cenvert old runlogs into new runlogs.
c
      implicit none

      real*8 
     & oblat,            ! observation latitude (deg).
     & oblon,            ! observation longitude (deg).
     & obalt,             ! observation altitude (km)
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
     & specname*57,      ! spectrum name
     & apf*2,            ! apodization function (e.g. BX N2, etc)
     & infile*24,
     & outfile*24 

      integer*4 ispec,nlhead,ncol,j,
     & lunw,
     & lunr,lunr2,      ! Logical unit number
     & istat,            ! status flag (0=success, 1=EOF)
     & iyr,              ! year
     & iset,             ! day of year
     & ifirst,           ! index of first spectral point in disk file
     & ilast,            ! index of last spectral point in disk file
     & possp,            ! Length of attached header in bytes
     & bytepw            ! Bytes per data word, usually 2 (I*2) or 4 (R*4)

      parameter (lunr=14,lunr2=15,lunw=16)

      col1=' '
      outfile='convert_runlog.out'
      infile='cruise94.grl'
      write(*,'(a)')'Name of old runlog:'
      read (*,'(a)') infile
      open(lunr,file=infile,status='old')
      read(lunr,*)nlhead,ncol
      do j=2,nlhead
         read(lunr,*)
      enddo

      open(lunw,file=outfile,status='unknown')

      write(lunw,*) ' 3   36'
      write(lunw,*) 'convert_runlog'
      write(lunw,'(a)')' Spectrum_File_Name                     Year'//
     &  '  Day   Hour   oblat    oblon   obalt    ASZA    POFF   AZIM'//
     &  '  OSDS    OPD   FOVI  FOVO  AMAL   IFIRST    ILAST'//
     &  '    DELTA_NU    POINTER  BPW ZOFF SNR  APF tins  pins'//
     &  '  hins   tout   pout  hout  sia    fvsi   wspd  wdir'//
     &  '  lasf    wavtkr  aipl'

      do ispec=1,9999999

          call read_runlog(lunr,col1,specname,iyr,iset,zpdtim,
     &    oblat,oblon,obalt,asza,zenoff,azim,osds,opd,fovi,fovo,amal,
     &    ifirst,ilast,graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     &    tout,pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)

          if (istat.ne.0) exit

          call write_runlog(lunw,col1,specname,iyr,iset,zpdtim,
     &    oblat,oblon,obalt,asza,zenoff,azim,osds,opd,fovi,fovo,amal,
     &    ifirst,ilast,graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     &    tout,pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)
      end do  !  ispec=1,9999999

      close(lunr)
      close(lunw)
      stop
      end
