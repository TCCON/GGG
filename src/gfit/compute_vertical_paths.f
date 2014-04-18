         subroutine compute_vertical_paths(ncell,zmin,z,d,vpath,nlev)
c
c  computes the *effective* slant path distances vertically above ZMIN
c  subject to the constraint that  Vertical_Column = SUM[VPATH(i).D(i)]
c
c  Inputs:
c     zmin         R*4  Minimum altitude (km)
c     nlev         I*4  Number of atmospheric levels
c     ncell        I*4  Number of non-atmospheric levels (cells)
c     z(nlev)      R*4  Level Altitudes (km)
c     d(nlev)      R*4  Level Densities (molecules/cm-3)
c
c  Output:
c     vpath(nlev)  R*4  Vector of effective vertical path distances
c
c  Notes:
c  Density is assumed to decrease exponentially with increasing altitude
c     d(z) = d(i-1) * exp[-((z-z(i-1))/(z(i)-z(i-1))*loge(d(i-1)/d(i))]
c  Note that d(z) = d(i-1) when z=z(i-1)  and  d(z) = d(i) when z=z(i)
c
c  Due to this exponential dependence,  the average density in the
c  vicinity of z(i) is greater than d(i).  To compensate for this,
c  the effective vertical paths must be greater than [Z(i+1)-Z(i-1)]/2
c  in order for the total vertical column to equal  VPATH.D
c
c  VPATH(i) = Integral [ (z-z(i-1)).d(z) ]  z=z(i-1),z(i)
c           + Integral [ (z(i+1)-z).d(z) ]  z=z(i),z(i+1)
c
c  VPATH(i) = [z(i)-z(i-1)]/loge(d(i-1)/d(i)).[1-d(i)/d(i-1)]
c           + [z(i+1)-z(i)]/loge(d(i)/d(i+1)).[1-d(i+1)/d(i)]

      implicit none
      integer*4 ncell,nlev,klev,jlev
      real*4 zmin,z(nlev),d(nlev),vpath(nlev),
     & logrp,dz,xo,xl

      call vmov(0.0,0,vpath,1,nlev)
c      write(*,*)'compute_columns: nlev=',nlev
      do klev=ncell+2,nlev
         if(z(klev).gt.zmin) go to 777
      end do
      klev=nlev
c      write(6,*) 'Warning: zmin exceeds z(nlev)',zmin,z(klev)
 777  dz=z(klev)-z(klev-1)
      xo=(zmin-z(klev-1))/dz
      if(d(klev).le.0.0) then
         logrp=0.0
      else
         logrp=log(d(klev-1)/d(klev))
      endif
c      write(6,*)'xo=',klev,xo,zmin,z(klev),z(klev-1),d(klev-1),d(klev)
      xl=logrp*(1.-xo)
      vpath(klev-1)=dz*(1.-xo)*
     &(1-xo-xl*(1+2*xo)/3+xl**2*(1+3*xo)/12+xl**3*(1+4*xo)/60)/2
      vpath(klev)  =dz*(1.-xo)*
     &(1+xo+xl*(1+2*xo)/3+xl**2*(1+3*xo)/12-xl**3*(1+4*xo)/60)/2
      do jlev=klev+1,nlev
        dz=z(jlev)-z(jlev-1)
        if(d(jlev).le.0.0) stop 'Non-positive density'
        logrp=log(d(jlev-1)/d(jlev))
        vpath(jlev-1)=vpath(jlev-1)+
     &  dz*(1.-logrp/3+logrp**2/12-logrp**3/60)/2
        vpath(jlev)  = dz*(1.+logrp/3+logrp**2/12+logrp**3/60)/2
      end do
      return
      end
