      subroutine setvmr(vmr,mgas,z,h2ovmr,co2vmr,nlev,iyr,iday,
     & uttime,oblat,oblon,obalt,tout,pout,hout,ztrop,zpbl)
c
c  Purpose:    Allows a vmr profile to be modified
c
c  Inputs:     Everything
c
c  Outputs:    vmr
c
      implicit none

      integer*4
     & mgas,       ! Maximum allowed number of gases
     & nlev, ilev, ! Number of atmospheric levels (70?)
     & iyr,        ! year (e.g. 2005)
     & iday        ! Day of year (1-366)

      real*4 
     & vmr(mgas,nlev),   ! buffer for vmr's
c     & p(nlev),          ! pressures of levels (atm.)
c     & t(nlev),          ! temperatures of levels (K)
     & z(nlev),          ! altitudes of levels (km)
     & h2ovmr(nlev),     ! H2O vmr profile from .mod file
     & co2vmr(nlev),     ! CO2 vmr profile from simulate_co2_vmr
     & co_vmr_max        ! maximum allowed CO vmr

      real*8 
     & uttime,           ! UT Time (hours)
     & oblat,            ! Observation Latitude (deg)
     & oblon,            ! Observation Longitude (deg)
     & obalt,            ! Observation Altitude (km)
     & tout,             ! Outside temperature (c)
     & pout,             ! Outside pressure (mbar)
     & hout,             ! Outside RH (%)
     & ztrop,            ! Tropopause  Altitude (km)
     & zpbl              ! PBL Altitude (km)

      if(1.eq.2) write(*,*)uttime,hout,oblon,obalt,tout,z  ! Avoids compiler warning for unused variables

c  Overwrite CO2 profile with simulation
c      write(*,*) 'ztrop=',ztrop
      if(ztrop.gt.0) then  ! fudge CO2 only if ztrop is non-zero.
          call simulate_co2_vmr(iyr,iday,uttime,oblat,oblon,
     &    ztrop,zpbl,z,co2vmr,nlev)
          do ilev=1,nlev
             vmr(2,ilev)=co2vmr(ilev)
          end do
      endif  ! if(ztrop.gt.0.0)

c  For ground-based observations, over-write apriori H2O vmr with NCEP model.
c  For HDO, the factor 0.15*(8.0+log10(h2ovmr(ilev))) adjusts for the effects
c  of isotopic fractionation.
c  For H2O = 0.03,  HDO=0.975*H2O (sea-level in tropics)
c  For H2O = 0.01,  HDO=0.900*H2O (sea-level at mid-latitude)
c  For H2O = 3E-06, HDO=0.375*H2O (tropopause)
      if(h2ovmr(1).gt.0.0) then ! An H2O profile was found in the .mod file
         if ( pout.gt.600.0 ) then ! ground-based observations
c            write(*,*) 'Replacing H2O & HDO vmrs with NCEP profiles'
            do ilev=1,nlev
              vmr(1,ilev)=h2ovmr(ilev)  ! H2O
              vmr(49,ilev)=h2ovmr(ilev)*0.16*(8.0+log10(h2ovmr(ilev))) ! HDO
            end do
         endif
      endif

c Limit tropospheric CO vmr to ~50 ppb in Southern hemisphere
c to 90 ppb in the tropics, and to 130 ppb in the northern hemisphere.
      co_vmr_max = 1.0E-09*(90 + 40*oblat/sqrt(200.0+oblat**2))
      do ilev=1,nlev
         if(z(ilev) .ge. ztrop) exit
         if(vmr(5,ilev).gt.co_vmr_max) vmr(5,ilev)=co_vmr_max
      end do

      return
      end
