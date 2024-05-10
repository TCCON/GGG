      subroutine apply_secular_trends(nlev,z,mgas,ngas,strend,vmrin,
     & ztrop_mod, lat_ref, alat_obs, date_obs,date_vmr,vmrout)
c
c  Purpose: 
c      Modifies the a priori profiles on a gas-by-gas basis
c      to account for the difference in time between the 
c      observation/model and the reference vmrs. This includes
c      the secular trend but not the seasonal cycle.
c
c  Inputs:
c      nlev               I*4   Number of levels
c      z(nlev)            R*4   Altitude (km) of levels
c      ngas               I*4   Number of gases
c      vmrin(mgas,nlev)   R*4   Reference VMR profile (from .vmr file)
c      ztrop_mod          R*8   Observation/Site tropopause altitude
c      lat_ref            R*8   Latitude of Reference vmrs
c      alat_obs           R*8   Latitude of site/observation
c      date_obs           R*8   Date of observation (.mod file)
c      date_vmr           R*8   Date of .vmr file)
c
c  Outputs:
c      vmrout(mgas,nlev)  R*4   Adjusted vmr profiles
c
c
c
      implicit none
      integer*4 ilev,nlev,mgas,ngas,jgas

      real*8 
     & date_obs,date_vmr,ztrop_mod,lat_ref,alat_obs,zobs,calc_aoa

      real*4 z(nlev),vmrin(mgas,nlev),vmrout(mgas,nlev),
     & strend(mgas),aoa,tdiff,tdmaoa

      lat_ref=lat_ref  ! avoid compiler warning

c  Apply secular trends
c      do jgas=1,mgas
c         strend(jgas)=0.0
c      end do
c      strend(2)=0.0052  ! CO2
c      strend(4)=0.001   ! N2O  (0.1% per year)
c      strend(5)=-.006   ! CO
c      strend(6)=0.003   ! CH4
c      strend(14)=-.01   ! HF
c      strend(42)=0.05   ! CHClF2 (HCFC-22)


c  Apply secular trends
      tdiff=sngl(date_obs-date_vmr)
      do ilev=1,nlev
         zobs=dble(z(ilev))
c         write(*,*)ilev,alat_obs,zobs,ztrop_mod
         aoa=sngl(calc_aoa(alat_obs,zobs,ztrop_mod))
         tdmaoa=tdiff-aoa
c         write(*,*)ilev,tdiff,aoa,tdmaoa,vmrin(2,ilev)
         do jgas=1,ngas
            vmrout(jgas,ilev)=vmrin(jgas,ilev)*
     &      (1+strend(jgas)*tdmaoa)
            if(jgas.eq.2) vmrout(jgas,ilev)=vmrout(jgas,ilev)*
     &      (1+(tdmaoa/155.0)**2)
            if(jgas.eq.6) vmrout(jgas,ilev)=vmrout(jgas,ilev)*
     &      (1.004-0.024*(tdmaoa+2.5)/sqrt(25.+(tdmaoa+2.5)**2))
            if(jgas.eq.14) vmrout(jgas,ilev)=vmrout(jgas,ilev)/
     &      (1.0+exp((-tdmaoa-16)/5.0))
            if(jgas.eq.50) vmrout(jgas,ilev)=1.5*vmrout(jgas,ilev)/
     &      (1.0+exp((-tdmaoa-4)/9.0))
         end do    !  do jgas=1,ngas
      end do     ! do ilev=1,nlev
c      write(74,*) date_obs,alat_obs,vmrout(2,16)

      return
      end
