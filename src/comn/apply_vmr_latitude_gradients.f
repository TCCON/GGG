      subroutine apply_vmr_latitude_gradients(nlev,z,mgas,ngas,gradlat,
     & refvmr,ztrop_mod,lat_ref,lat_obs,aprvmr)
c
c  Purpose:  Modifies the a priori profiles on a gas-by-gas basis to
c            account for the difference in latitude between the site
c            and the reference latitude for the a priori profile (35N).
c
c  Inputs:
c          nlev               I*4   Number of levels
c          z(nlev)            R*4   Altitude (km) of levels
c          mgas               I*4   Number of gases
c          ngas               I*4   Number of gases
c          gradlat(mgas)      R*4   Latitude gradients
c          refvmr(mgas,nlev)  R*4   Reference VMR profile (from .vmr file)
c          ztrop_mod          R*8   Observation/Site tropopause altitude
c          lat_ref            R*8   Latitude of Reference vmrs
c          lat_obs            R*8   Latitude of site/observation
c
c  Outputs:
c          apvmr(mgas,nlev)   R*4   Adjusted vmr profiles
c
c
c  Comments:
c    
c   In the upper stratosphere, gas distributions are assumed symmetrical about equator,
c   At the surface, gas distributions are assumed anti-symmetric about equator.
c   At intermediate altitudes the profiles are interpolated between these limiting
c   behaviors.

c  Assume an initial vmr profile vmr(z,igas) appropriate for lat_ref
c  Assume that at the surface the vmr varies with latitude as
c     vmr(lat) = vmr(0)*(1+gradlat*(lat/15)/sqrt(1+(lat/15)**2))
c  This function varies smoothly from:
c     vmr(0)*(1+gradlat) at the NP
c     vmr(0)*(1+gradlat/2) at the 15N
c     vmr(0) at the  equator
c     vmr(0)*(1-gradlat/2) at the 15S
c     vmr(0)*(1-gradlat) at the SP
c  with the maximum change at the equator.
c
c     vmr(lat_ref) = vmr(0)*(1+gradlat*(lat_ref/15)/dsqrt(1+(lat_ref/15)**2))
c     vmr(lat_obs) = vmr(0)*(1+gradlat*(lat_obs/15)/dsqrt(1+(lat_obs/15)**2))
c
      implicit none
      integer*4 ilev,nlev,mgas,ngas,jgas

      real*8 ztrop_mod,lat_ref,lat_obs,xobs,xref,fr

      real*4 z(nlev),gradlat(mgas),refvmr(mgas,nlev),aprvmr(mgas,nlev)

c  Inter-hemispheric gradients
c  Positive values imply more in the NH than the SH than
c  would be expected based solely on the age of the air.
c       do jgas = 1,mgas
c         gradlat(jgas)=0.0   ! default
c       end do
c       gradlat(2) =0.0000 ! CO2
c       gradlat(3) =0.20   ! O3
c       gradlat(5) =0.28   ! CO
c       gradlat(6) =0.033  ! CH4
c       gradlat(10)=0.250  ! NO2
c       gradlat(11)=0.200  ! NH3
c       gradlat(12)=0.100  ! HNO3
c       gradlat(20)=0.200  ! H2CO
c       gradlat(28)=0.100  ! HCN
c       gradlat(29)=0.200  ! CH3F
c       gradlat(30)=0.200  ! CH3Cl
c       gradlat(31)=0.200  ! CF4
c       gradlat(32)=0.200  ! CCl2F2
c       gradlat(33)=0.200  ! CCl3F
c       gradlat(34)=0.200  ! CH3CCl3
c       gradlat(35)=0.200  ! CCl4
c       gradlat(38)=0.300  ! C2H6
c       gradlat(39)=0.300  ! C2H4
c       gradlat(40)=0.300  ! C2H2
c       gradlat(42)=0.200  ! CHClF2
c       gradlat(44)=0.200  ! CH3Br
c       gradlat(46)=0.200  ! HCOOH
c       gradlat(48)=0.200  ! CHCl2F
c       gradlat(50)=0.300  ! SF6
c       gradlat(51)=0.300  ! F113
c       gradlat(53)=0.200  ! F142b
c       gradlat(56)=0.200  ! CH3OH
c       gradlat(58)=0.200  ! CH3CHO
c       gradlat(59)=0.200  ! CH3CN
c       gradlat(61)=0.300  ! NF3
c       gradlat(65)=0.200  ! CHF3
c       gradlat(66)=0.200  ! f141b
c       gradlat(67)=0.200  ! CH3COOH
c       gradlat(70)=0.500  ! C3H8   

      do jgas=1,ngas
         xref=gradlat(jgas)*(lat_ref/15)/dsqrt(1+(lat_ref/15)**2)
         xobs=gradlat(jgas)*(lat_obs/15)/dsqrt(1+(lat_obs/15)**2)
         do ilev=1,nlev
            fr=1.0/(1.0+(z(ilev)/ztrop_mod)**2)
            aprvmr(jgas,ilev)=refvmr(jgas,ilev)*
     &      sngl((1+fr*xobs)/(1+fr*xref))
c            if(jgas.eq.70) write(*,'(4f11.4,2e12.4)')xref,xobs,
c     &      z(ilev),fr,refvmr(jgas,ilev),aprvmr(jgas,ilev)
         end do
      end do    !  do jgas=1,ngas

      return
      end
