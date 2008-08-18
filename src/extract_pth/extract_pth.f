c  Program extract_pth
c  Interpolates into a NCEP model files to determine
c  the T/P/H at the altitude of the site, and then output these
c  values along with the corresponding values from the runlog.
c  Useful for comparing the weather station data with NCEP.

      implicit none
      integer*4 i,iwas,lunm,lunw,ncol,nlhead
      integer*4
     & lunr,             ! Logical unit number
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
     & snr,              ! Signal-to-Noise Ratio (dimensionless)
     & tins,             ! Inside temperature
     & pins,             ! Inside pressure
     & hins,             ! Inside humidity
     & tout,             ! Outside temperature
     & pout,             ! Outside pressure
     & hout,             ! Outside humidity
     & sia,              ! Solar Intensity (Average)
     & sis,              ! Solar Intensity (SD)
     & aipl,             ! Airmass-Independent Path Length (km)
     & lasf,             ! Laser Frequency (e.g. 15798 cm-1)
     & wavtkr           ! suntracker frequency (active tracking)

      character
     & col1*1,           ! first column of runlog record
     & runlab*21,        ! spectrum name
     & apf*2             ! apodization function (e.g. BX N2, etc)

      parameter (lunr=14,lunm=15,lunw=16)
      character modfile*14, modwas*14

      real*4 roc,ecc2,alat,gs,zz,pfact,ptrop,pmod,tmod,hmod,
     &   scht,log1pxox,
     &   p0,t0,z0,m0,h0,p1,t1,z1,m1,h1

      iwas=0
      modwas='            '
      open(lunr, file='/home/toon/ggg/runlogs/gnd/paIn_1.0fg.grl',
     & status='old')
      read(lunr,*)
      open(lunw,file='extract_pth.out', status='unknown')
      write(lunw,*)2,9
      write(lunw,*)' yy doy  hh  tout tmod  pout pmod  hout hmod'
      
      do i=1,999999
       call read_runlog(lunr,col1,runlab,iyr,iset,zpdtim,
     & oblat,oblon,obalt,asza,zenoff,opd,fovi,fovo,amal,ifirst,
     & ilast,graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     & tout,pout,hout,lasf,wavtkr,sia,sis,aipl,istat)
      if(istat.ne.0) exit
      modfile=runlab(:10)//'.mod'
      if(modfile.ne.modwas) then
         open(unit=lunm,file='/home/toon/ggg/models/gnd/'//modfile,
     &   status='old')
         read(lunm,*)nlhead,ncol
         read(lunm,*)roc,ecc2,alat,gs,zz,pfact,ptrop
         read(lunm,*)
         read(lunm,*)
         read(lunm,*)p0,t0,z0,m0,h0
         read(lunm,*)p1,t1,z1,m1,h1
         close(lunm)
         write(*,*)modwas, i-iwas
         iwas=i
         modwas=modfile
      endif
      tmod=t0 + (t1-t0)*(obalt-z0)/(z1-z0)
      hmod=h0 + (h1-h0)*(obalt-z0)/(z1-z0)
      pmod=p0*exp(-(obalt-z0)/8)
      scht=8.3144*t0/m0/gs/log1pxox(tmod/t0-1)
      pmod=p0*exp(-(obalt-z0)/scht)
      write(lunw,'(2i4,7f9.3)') iyr,iset,zpdtim,
     & tout+273.16,tmod,pout,pmod,hout,hmod
      end do
      write(*,*)i-1
      close(lunr)
      stop
      end
