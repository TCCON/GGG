      subroutine apply_vmr_latitude_gradients(nlev,z,mgas,refvmr,
     & ztrop_mod, lat_ref, lat_obs, aprvmr)
c
c  Purpose:  Modifies the a priori profiles on a gas-by-gas basis to
c            account for the difference in latitude between the site
c            and the reference latitude for the a priori profile (35N).
c
c  Inputs:
c          nlev               I*4   Number of levels
c          z(nlev)            R*4   Altitude (km) of levels
c          ngas               I*4   Number of gases
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
c   In the middle stratosphere, gas distributions are assumed symmetrical about equator,
c   At the surface, gas distributions are assumed anti-symmetric about equator.
c   At intermediate altitudes the profiles are interpolated between these limiting
c   behaviors.

c  Assume an initial vmr profile vmr(z,igas) appropriate for lat_ref
c  Assume that at the surface the vmr varies with latitude as
c     vmr(lat) = vmr(0)*(1+gradlat*(lat/15)/sqrt(1+(lat/15)**2))
c  This function varies smoothly from:
c     vmr(0)*(1+gradlat) at the NP
c     vmr(0) at the  equator
c     vmr(0)*(1-gradlat) at the SP
c  with the maximum change at the equator.
c
c     vmr(lat_ref) = vmr(0)*(1+gradlat*(lat_ref/15)/sqrt(1+(lat_ref/15)**2))
c     vmr(lat_obs) = vmr(0)*(1+gradlat*(lat_obs/15)/sqrt(1+(lat_obs/15)**2))
c
      implicit none
      integer*4 ilev,nlev,mgas,kgas,jgas

      real*8 gradlat(6),ztrop_mod,lat_ref,lat_obs,xobs,xref,fr

      real*4 z(nlev),refvmr(mgas,nlev),aprvmr(mgas,nlev)

c  Inter-hemispheric gradients
      gradlat(1)=0.00    ! H2O
      gradlat(2)=0.0015   ! CO2
      gradlat(3)=0.20    ! O3
      gradlat(4)=0.000   ! N2O
      gradlat(5)=0.25    ! CO
      gradlat(6)=0.035   ! CH4
      kgas=6

      do jgas=1,kgas
         xref=gradlat(jgas)*(lat_ref/15)/sqrt(1+(lat_ref/15)**2)
         xobs=gradlat(jgas)*(lat_obs/15)/sqrt(1+(lat_obs/15)**2)
         do ilev=1,nlev
c           if(z(ilev).lt.ztrop_mod) then
c             fr=1.0-(z(ilev)/ztrop_mod)**2
             fr=1.0/(1.0+(z(ilev)/ztrop_mod)**2)
             aprvmr(jgas,ilev)=refvmr(jgas,ilev)*(1+fr*xobs)/(1+fr*xref)
c           else
c             aprvmr(jgas,ilev)=refvmr(jgas,ilev)
c           endif
         end do
      end do    !  do jgas=1,kgas
      return
      end
