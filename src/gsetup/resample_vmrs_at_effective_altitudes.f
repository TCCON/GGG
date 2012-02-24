      subroutine resample_vmrs_at_effective_altitudes(nlev,z,mgas,ngas,
     & refvmr,ztrop_mod,ztrop_vmr,obslat_mod,reflat_vmr,adjvmr)
c  Reads and interpolates (possibly with a vertical stretch) an initial vmr set
c  (SETNAME) onto the vertical grid defined by Z(NLAY) and places the result
c  into array VMR.

c  Inputs:
c     nlev              I*4    Number of levels
c     z(nlev)           R*4    Altitudes of levels (km)
c     mgas              I*4    Declared first dimension of VMR matrix
c     ngas              I*4    Actual number of gases
c     refvmr(mgas,nlev) R*4    Gas VMRs
c     ztrop_mod         R*8    Tropopause altitude from .mod file
c     ztrop_vmr         R*8    Tropopause altitude from .vmr file
c     obslat_mod        R*8    Observation/Site Latitude from .mod file
c     reflat_vmr        R*8    Latitude of reference vmr set
c
c Outputs:
c     adjvmr(mgas,nlev) R*4    Gas VMRs resampled in altitude


      implicit none

      INTEGER*4 klev,nlev,jgas,mgas,ngas,jj
      real*4 z(nlev),refvmr(mgas,nlev),adjvmr(mgas,nlev),fr
      real*8 zeff,ztrop_mod,ztrop_vmr,stuclat,obslat_mod,reflat_vmr
c==================================================================
      do klev=1,nlev
c         write(*,*)klev,nlev,z(klev),ztrop_mod,ztrop_vmr
         if(z(klev).lt.ztrop_mod) then
            zeff=z(klev)*ztrop_vmr/ztrop_mod  ! troposphere
         else
            zeff=z(klev)+(ztrop_vmr-ztrop_mod)*
     &      exp(-(z(klev)-ztrop_mod)/10.0)       ! stratosphere
            zeff=zeff-exp(-((obslat_mod-stuclat)/25.0)**4)
     &     *3.5*ztrop_mod*(z(klev)/ztrop_mod-1)**2
     &     *exp(-(z(klev)-ztrop_mod)/9.0)
         endif
c         write(*,*) 'zeff=',zeff
         if(zeff.gt.z(nlev)) zeff=z(nlev)

c   Try to find jj such that   z(jj-1) < zeff < z(jj)
         do jj=2,nlev-1
            if(z(jj).gt.zeff) exit
         end do
         fr=(zeff-z(jj-1))/(z(jj)-z(jj-1))
         do jgas=1,ngas
         adjvmr(jgas,klev)=fr*refvmr(jgas,jj)+(1.-fr)*refvmr(jgas,jj-1)
         end do
      end do  ! klev=1,nlev

      return
      end
