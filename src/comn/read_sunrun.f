      subroutine read_sunrun(lun,col1,specname,object,tcorr,oblat,
     &   oblon,obalt,tins,pins,hins,tout,pout,hout,sia,fvsi,
     &   wspd,wdir,fmin,fmax,fsf,lasf,wavtkr,aipl,tel_mag,istat)

c  Reads a single record from a sunrun file, which must already
c  have been opened to logical unit number LUN.  The goal of this
c  subroutine is to hide all the of the code that depends on the
c  sunrun format into pair of subroutines (writsunrun & readsunrun).
c
c  This has two advantages:
c  1) It simplifies the calling programs
c  2) It means that if the sunrun format is changed
c  in the future, only writsunrun & readsunrun need be changed,
c  not the dozen main programs that read the sunrun.
c
c  Input:
c    lun   Logical Unit number of file to be read
c    everything else

c  Outputs:
c    istat

      implicit none

      integer*4
     & lun,              ! Logical unit number
     & object,           ! Heavenly object (Moon=1, Sun=2)
     & istat             ! status flag (0=success, 1=EOF)

      real*8
     & tcorr,            ! Time Correction (s)
     & oblat,            ! Observation latitude (deg).
     & oblon,            ! Observation longitude (deg).
     & obalt,            ! Observation altitude (km)
     & tins,             ! Inside temperature (C)
     & pins,             ! Inside pressure (mbar)
     & hins,             ! Inside humidity (%)
     & tout,             ! Outside temperature (C)
     & pout,             ! Outside pressure (mbar)
     & hout,             ! Outside humidity (%)
     & sia,              ! Solar Intensity Average (arbitrary units)
     & fvsi,             ! Fractional Variation in Solar Intensity
     & wspd,             ! Wind speed (m/s)
     & wdir,             ! Wind Direction (deg)
     & aipl,             ! Airmass-Independent Path Length (km)
     & tel_mag,          ! Telescope Magnification (dimensionless)
     & fmin,             ! Start frequency (cm-1)
     & fmax,             ! End frequency (cm-1)
     & fsf,              ! Frequency Stretch Factor
     & lasf,             ! Laser Frequency (e.g. 15798 cm-1)
     & wavtkr            ! Suntracker frequency (active tracking)

      character
     & col1*1,           ! first column
     & specname*(*)      ! spectrum name

1     read(lun,34,end=99) col1,specname,object,tcorr,oblat,oblon,
     & obalt,tins,pins,hins,tout,pout,hout,sia,fvsi,wspd,wdir,
     & fmin,fmax,fsf,lasf,wavtkr,aipl,tel_mag 
      if( col1.eq.':' .or. col1.eq.';' ) go to 1
c      write(*,*) specname

 34   format(a1,a57,1x,i2,f8.0,f9.4,f10.4,
     & f7.3,f6.1,f8.2,f6.1,f6.1,f8.2,f6.1,f7.1,f7.4,f6.1,f6.0,
     & 1x,2f7.0,f11.8,f11.3,f7.0,f7.3,f6.2)
      istat=0
      return
99    istat=1
      return
      end
