      subroutine setvmr(vmr,mgas,z,t,p,h2ovmr,nlev,iyr,iday,
     & uttime,oblat,oblon,obalt,tout,pout,hout,ptrop_atm,ppbl_atm)
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
     & p(nlev),          ! pressures of levels (atm.)
     & t(nlev),          ! temperatures of levels (K)
     & z(nlev),          ! altitudes of levels (km)
     & co2vmr(nlev),     ! H2O vmr profile from .mod file
     & h2ovmr(nlev)      ! H2O vmr profile from .mod file

      real*8 
     & uttime,           ! UT Time (hours)
     & oblat,            ! Observation Latitude (deg)
     & oblon,            ! Observation Longitude (deg)
     & obalt,            ! Observation Altitude (km)
     & tout,             ! Outside temperature (c)
     & pout,             ! Outside pressure (mbar)
     & hout,             ! Outside RH (%)
     & ptrop_atm,        ! Tropopause pressure (atm)
     & ppbl_atm          ! PBL pressure (atm)

      if(1.eq.2) write(*,*)uttime,hout,oblon,obalt,tout,t,z  ! Avoids compiler warning for unused variables

c  Overwrite CO2 profile with simulation
      if(ptrop_atm.gt.0) then  ! fudge CO2 only if ptrop is non-zero.
          call simulate_co2_vmr(iyr,iday,uttime,oblat,oblon,
     &    ptrop_atm,ppbl_atm,z,p,co2vmr,nlev)
          do ilev=1,nlev
             vmr(2,ilev)=co2vmr(ilev)
          end do
      endif  ! if(ptrop_atm.gt.0.0)


c  Substitute NCEP model H2O profile for ground-based observations.
c  The factor 0.16*(8.0+log10(h2ovmr(ilev))) represents
c  the isotopic fractionation of HDO/H2O.
      if(h2ovmr(1).gt.0.0) then ! An H2O profile was found in the .mod file
         if ( pout.gt.600.0 ) then ! ground-based observations
            write(*,*) 'Replacing H2O & HDO vmrs with NCEP profiles'
            do ilev=1,nlev
              vmr(1,ilev)=h2ovmr(ilev)  ! H2O
              vmr(49,ilev)=h2ovmr(ilev)*0.16*(8.0+log10(h2ovmr(ilev))) ! HDO
            end do
         endif
      endif

      return
      end
