      subroutine write_sunrun(lun,col1,specname,obj,tcorr,oblat,
     &   oblon,obalt,tins,pins,hins,tout,pout,hout,nus,nue,fsf,
     &   lasf,wavtkr,sia,sis,aipl,tel_mag,istat)
c
c  Writes a single record to a sunrun file.
c  The file must already have been opened.
c  The goal of this subroutine is to hide all the of the code
c  that depends on the sunrun format into pair of subroutines
c  (writsunrun & readsunrun).
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
     & obj,              ! Heavenly object (Moon=1, Sun=2)
     & istat             ! status flag (0=success, 1=EOF)

      real*8
     & tcorr,            ! Time correction (s)
     & oblat,            ! Observation latitude (deg).
     & oblon,            ! Observation longitude (deg).
     & obalt,            ! Observation altitude (km)
     & tins,             ! Inside temperature (C)
     & pins,             ! Inside pressure (mbar)
     & hins,             ! Inside humidity (%)
     & tout,             ! Outside temperature (C)
     & pout,             ! Outside pressure (mbar)
     & hout,             ! Outside humidity (%)
     & sia,              ! Solar Intensity (Average)
     & aipl,             ! Airmass-Independent Path Length (km)
     & tel_mag,          ! Telescope Magification (dimensionless)
     & sis,              ! Solar Intensity (SD)
     & nus,              ! Start frequency (cm-1)
     & nue,              ! End frequency (cm-1)
     & fsf,              ! Frequency Stretch Factor
     & lasf,             ! Laser Frequency (e.g. 15798 cm-1)
     & wavtkr            ! Suntracker frequency (active tracking)

      character
     & col1*1,           ! first column
     & specname*21         ! spectrum name

      write(lun,34,err=99) col1,specname,obj,tcorr,oblat,oblon,obalt,
     & tins,pins,hins,tout,pout,hout,nus,nue,fsf,lasf,wavtkr,sia,sis,
     & aipl,tel_mag
 34   format(a1,a21,1x,i2,f8.0,f9.4,f10.4,f7.3,f6.1,f8.2,f6.1,f6.1,f8.2,
     & f6.1,1x,2f7.0,f11.8,f11.3,f7.0,2f6.1,f7.3,f6.2)
      istat=0
      return
99    istat=1
      return
      end
