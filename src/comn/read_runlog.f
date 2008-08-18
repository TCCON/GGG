      subroutine read_runlog(lun,col1,runlab,iyr,iset,zpdtim,
     & oblat,oblon,obalt,asza,zenoff,opd,fovi,fovo,amal,ifirst,
     & ilast,graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     & tout,pout,hout,lasf,wavtkr,sia,sis,aipl,istat)
c
c  Reads a single record from the runlog file.
c  File must already have been opened.
c  The purpose of this subroutine is to hide all the
c  of the code that depends on the runlog format into
c  a single subroutine. This has two advantages:
c  1) It simplifies the calling programs
c  2) It means that if the runlog format is changed
c  in the future, only read_runlog subroutine needs to be changed,
c  not the dozen main programs that read the runlog.
c
c
c  Input:
c    lun   Logical Unit number of file to be read

c  Outputs:
c    everything else

      implicit none
      integer*4 lnbc,lr,
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
     & wavtkr,           ! suntracker frequency (active tracking)
     & zerr,             ! Extra parameter in ATMOS tab-delimited runlogs
     & scalf             ! Extra parameter in ATMOS tab-delimited runlogs

      character
     & col1*1,           ! first column of runlog record
     & record*400,       ! runlog record
     & runlab*(*),       ! spectrum name
     & apf*2             ! apodization function (e.g. BX N2, etc)

1      read(lun,'(a1,a)',end=99) col1,record
      if( col1.eq.':') go to 1
      if(col1.ne.'-' .and. col1.ne.'+' .and. col1.ne.' ') then
         record=col1//record   ! Runlog is the old format
      endif
      lr=lnbc(record)
c      write(*,*)'read_runlog: lr= ', lr
c
c Note: ASCI character 9 is a horizontal tab.
      if(index(record,char(9)).gt.0) then    ! TAB delimited (e.g. ATMOS)
        read(record,*) runlab,iyr,iset,zpdtim,oblat,oblon,
     &  obalt,asza,zenoff,opd,fovi,fovo,amal,ifirst,ilast,graw,
     &  possp,bytepw,zoff,zerr,snr,scalf,apf
      elseif(lr.le.233) then     ! SPACE delimited (e.g. MkIV)
        read(record,332,err=97) runlab,iyr,iset,zpdtim,oblat,oblon,
     &  obalt,asza,zenoff,opd,fovi,fovo,amal,ifirst,ilast,graw,possp,
     &  bytepw,zoff,snr,apf,tins,pins,hins,tout,pout,hout,lasf,wavtkr,
     &  sia,sis,aipl
 332  format(a21,1x,2i4,f8.4,f8.3,f9.3,2f8.3,f7.0,f7.2,3f6.0,2i8,f15.11,
     & i8,i3,1x,2f5.0,1x,a2,2(f6.0,f8.0,f5.0),f10.0,f7.0,2f6.1,f7.3)
      else                 ! New GDS-format runlog
        read(record,333,err=98) runlab,iyr,iset,zpdtim,oblat,oblon,
     &  obalt,asza,zenoff,opd,fovi,fovo,amal,ifirst,ilast,graw,possp,
     &  bytepw,zoff,snr,apf,tins,pins,hins,tout,pout,hout,lasf,wavtkr,
     &  sia,sis,aipl
 333  format(a34,2x,2i4,f8.4,f8.3,f9.3,2f8.3,f7.0,f7.2,3f6.0,2i8,f15.11,
     & i8,i3,1x,2f5.0,1x,a2,2(f6.0,f8.0,f5.0),f10.0,f7.0,2f6.1,f7.3)
      endif
      istat=0
      return
97    istat=1
      return
98    istat=2
      return
99    istat=3
      return
      end
