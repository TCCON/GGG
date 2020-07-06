      subroutine resample_vmrs_at_effective_altitudes(nlev,z,mgas,ngas,
     & itcz_lat, itcz_width,
     & refvmr,ztrop_mod,ztrop_vmr,obslat_mod,reflat_vmr,adjvmr)
c  Reads and interpolates (possibly with a vertical stretch) an initial vmr set
c  (SETNAME) onto the vertical grid defined by Z(NLAY) and places the result
c  into array VMR.

c  Inputs:
c     nlev              I*4    Number of levels
c     z(nlev)           R*4    Altitudes of levels (km)
c     mgas              I*4    Declared first dimension of VMR matrix
c     ngas              I*4    Actual number of gases
c     itcz_lat          R*8    Center latitude of ITCZ (deg)
c     itcz_width        R*8    Width of ITCZ (deg)
c     refvmr(mgas,nlev) R*4    Gas VMRs
c     ztrop_mod         R*8    Tropopause altitude from .mod file
c     ztrop_vmr         R*8    Tropopause altitude from .vmr file
c     obslat_mod        R*8    Observation/Site Latitude from .mod file
c     reflat_vmr        R*8    Latitude of reference vmr set
c
c Outputs:
c     adjvmr(mgas,nlev) R*4    Gas VMRs resampled in altitude
c
c  In the troposphere    zeff=z(klev)*ztrop_vmr/ztrop_mod
c  So at the tropopause  z=ztrop_mod, so  zeff=ztrop_vmr

c  In the stratosphere  zeff=ztrop_vmr at the tropopause,
c  and heads asymptotically to  zeff = z(klev) at high altitudes
c
c  So the *difference* between zeff and z(ilev) is zero at
c  sea level, increases linearly up to the tropopause, then
c  decreases exponentially back to zero.
c
c  In the tropics the is an extra uplift (a reduction of zeff)
c  by a factor dz^2.exp(-dz.10)

      implicit none

      INTEGER*4 klev,nlev,jgas,mgas,ngas,jj
      real*4 z(nlev),refvmr(mgas,nlev),adjvmr(mgas,nlev),fr
      real*8 zeff,ztrop_mod,ztrop_vmr,dztr,
c     & exdztr,
     & itcz_lat,itcz_width,obslat_mod,reflat_vmr,rdum

      rdum=reflat_vmr  ! Prevent compiler warning (unused variable)

      do klev=1,nlev
c         write(*,*)klev,nlev,z(klev),ztrop_mod,ztrop_vmr
         if(z(klev).lt.ztrop_mod) then      ! troposphere
            zeff=z(klev)*ztrop_vmr/ztrop_mod
         else                               ! stratosphere
            dztr=z(klev)-ztrop_mod  ! distance above trop
            zeff=z(klev)+exp(-(z(klev)-ztrop_mod)/10.)*
     &      (ztrop_vmr-ztrop_mod-3.5*ztrop_mod*(z(klev)/ztrop_mod-1)**2*
     &      exp(-((obslat_mod-itcz_lat)/(itcz_width+10))**4))
         endif
c         write(*,*) 'zeff=',zeff
         if(zeff.gt.z(nlev)) zeff=z(nlev)

c   Try to find jj such that   z(jj-1) < zeff < z(jj)
         do jj=2,nlev-1
            if(z(jj).gt.zeff) exit
         end do
         fr=sngl(zeff-z(jj-1))/(z(jj)-z(jj-1))
         do jgas=1,ngas
            adjvmr(jgas,klev)=fr*refvmr(jgas,jj)+
     &      (1.-fr)*refvmr(jgas,jj-1)
         end do
      end do  ! klev=1,nlev

      return
      end
