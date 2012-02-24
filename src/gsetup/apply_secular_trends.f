      subroutine apply_secular_trends(nlev,z,mgas,vmrin,
     & ztrop_mod, lat_ref, alat_obs, date_obs,date_vmr,vmrout)
c
c  Purpose:  Modifies the a priori profiles on a gas-by-gas basis
c            to account for the difference in time between the 
c            observation/model and the reference vmrs. This includes
c            the secular trend but not the seasonal cycle.
c
c  Inputs:
c          nlev              I*4   Number of levels
c          z(nlev)           R*4   Altitude (km) of levels
c          ngas              I*4   Number of gases
c          vmrin(mgas,nlev)  R*4   Reference VMR profile (from .vmr file)
c          ztrop_mod         R*8   Observation/Site tropopause altitude
c          lat_ref           R*8   Latitude of Reference vmrs
c          alat_obs          R*8   Latitude of site/observation
c          date_obs,         R*8   Date of observation (.mod file)
c          date_vmr,         R*8   Date of .vmr file)
c
c  Outputs:
c          vmrout(mgas,nlev)   R*4   Adjusted vmr profiles
c
c
c
      implicit none
      integer*4 ilev,nlev,mgas,kgas,jgas
      parameter (kgas=6)

      real*8 trend(kgas),
     & date_obs,date_vmr,ztrop_mod,lat_ref,alat_obs,
     & tdiff,zobs

      real*4 z(nlev),vmrin(mgas,nlev),vmrout(mgas,nlev),
     & calc_aoa,aoa

      lat_ref=lat_ref  ! avoid compiler warning

c  Apply secular trends
      trend(1)=0.00    ! H2O
      trend(2)=0.0054   ! CO2
      trend(3)=0.00    ! O3
      trend(4)=0.001   ! N2O  (0.1% per year)
      trend(5)=-.005   ! CO
      trend(6)=0.004   ! CH4

      tdiff=date_obs-date_vmr
      do ilev=1,nlev
         zobs=z(ilev)
         aoa=calc_aoa(alat_obs,zobs,ztrop_mod)
         do jgas=1,kgas
             vmrout(jgas,ilev)=vmrin(jgas,ilev)*
     &                           (1+trend(jgas)*(tdiff-aoa))
         end do
      end do    !  do jgas=1,kgas

      return
      end
