      subroutine compute_vertical_paths(ncell,zmin,z,d,vpath,nlev)
c
c  Computes vector of *effective* path distances (VPATH) vertically above
c  ZMIN subject to the constraint:   Vertical_Column = SUM[VPATH(i).D(i)]
c
c  Inputs:
c     ncell        I*4  Number of non-atmospheric levels (cells)
c     zmin         R*4  Minimum altitude (km)
c     nlev         I*4  Number of atmospheric levels
c     z(nlev)      R*4  Level Altitudes (km)
c     d(nlev)      R*4  Level Densities (molecules/cm-3)
c
c  Output:
c     vpath(nlev)  R*4  Vector of effective vertical path distances
c
c  Simplest aproach would be to define the layer boundary altitudes
c         Z+ = (z(i+1)+z(i))/2 
c         Z- = (z(i)+z(i-1))/2 
c    vpath(i) = Z+ - Z- =  [z(i+1)-z(i-1)]/2
c  but this would assume a constant density d(i) between Z- and Z+,
c  ignoring the exponential decrease of density with altitude.
c
c  Due to this exponential dependence, the average density in the
c  vicinity of z(i) is greater than d(i).  To compensate for this,
c  the effective vertical paths must be greater than [Z(i+1)-Z(i-1)]/2
c  in order for the total vertical column to equal  SUM[VPATH(i).D(i)]
c
c  Assume a linear trapezoidal weighting of the density and path
c  length associated with each level. So at z=z(i) 100% of the
c  path is associated with VPATH(i). Above and below z(i), 
c  the contribution decreases linearly, reaching 0% at z(i+1)
c  and at z(i-1).
c
c  VPATH(i).d(i) = Integral [ (z-z(i-1)).d(z) ]  z=z(i-1),z(i)
c                + Integral [ (z(i+1)-z).d(z) ]  z=z(i),z(i+1)
c
c  For z(i-1) < z < z(i)
c      d(z) = d(i-1) * Exp[-(z-z(i-1))/(z(i)-z(i-1))*loge(d(i-1)/d(i))]
c      d(z) = d(i-1) * Exp[-(z-z(i-1))/h]   h = (z(i)-z(i-1))/loge(d(i-1)/d(i))
c  Note that d(z) = d(i-1) when z=z(i-1)  and  d(z) = d(i) when z=z(i)
c
c  For the slant column between z(i-1) and z(i), 
c     vpath(i-1) = vpath(i-1) + Integral d(z).(z(i)-z)/(z(i)-z(i-1)) dz
c     vpath(i)   = vpath(i)   + Integral d(z).(z-z(i-1))/(z(i)-z(i-1)) dz
c
c  Integral [(z0-z) Exp[-(z-z1)/h] = +h.(h+z-z0).Exp[-(z-z1)/h] z=z1,z0
c    =  h.h.Exp[-(z0-z1)/h] - h(h+z1-z0)

c  Integral [(z-z1) Exp[-(z-z1)/h] = -h.(h+z-z1).Exp[-(z-z1)/h] z=z1,z0
c    =  h.(z1-z0).Exp[-(z1-z0)/h]

c  VPATH(i).d(i) = [z(i)-z(i-1)]/loge(d(i-1)/d(i)).[1-d(i)/d(i-1)]
c                + [z(i+1)-z(i)]/loge(d(i)/d(i+1)).[1-d(i+1)/d(i)]

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
     & (1-xo-xl*(1+2*xo)/3+xl**2*(1+3*xo)/12+xl**3*(1+4*xo)/60)/2
      vpath(klev)  =dz*(1.-xo)*
     & (1+xo+xl*(1+2*xo)/3+xl**2*(1+3*xo)/12-xl**3*(1+4*xo)/60)/2
      do jlev=klev+1,nlev
         dz=z(jlev)-z(jlev-1)
         if(d(jlev).le.0.0) then
            write(*,*) jlev, z(jlev), d(jlev)
            stop 'compute_vertical_paths: Non-positive density'
         endif
         logrp=log(d(jlev-1)/d(jlev))
         vpath(jlev-1)=vpath(jlev-1)+
     &   dz*(1.-logrp/3+logrp**2/12-logrp**3/60)/2
         vpath(jlev) = dz*(1.+logrp/3+logrp**2/12+logrp**3/60)/2
      end do
      return
      end
