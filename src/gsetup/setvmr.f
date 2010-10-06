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
     & co2vmr(nlev)      ! CO2 vmr profile from simulate_co2_vmr

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
c  The factor 0.16*(8.0+log10(h2ovmr(ilev))) represents
c  the isotopic fractionation of HDO/H2O.
      if(h2ovmr(1).gt.0.0) then ! An H2O profile was found in the .mod file
         if ( pout.gt.600.0 ) then ! ground-based observations
c            write(*,*) 'Replacing H2O & HDO vmrs with NCEP profiles'
            do ilev=1,nlev
              vmr(1,ilev)=h2ovmr(ilev)  ! H2O
              vmr(49,ilev)=h2ovmr(ilev)*0.16*(8.0+log10(h2ovmr(ilev))) ! HDO
            end do
         endif
      endif

      return
      end
