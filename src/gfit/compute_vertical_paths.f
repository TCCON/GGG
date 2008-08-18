         subroutine compute_vertical_paths(zmin,z,d,vpath,nlev)
c
c  computes the slant path distances vertically above ZMIN
c
c  Inputs:
c     zmin     R*4  Minimum altitude (km)
c     z(nlev)  R*4  Level Altitudes
c     d(nlev)  R*4  level Densities
c
c  Output:
c     vpath(nlev)  R*4  Vector of vertical path distances

      implicit none
      integer*4 nlev,klev,jlev
      real*4 zmin,z(nlev),d(nlev),vpath(nlev),
     &logrp,dz,xo,xl

      call vmov(0.0,0,vpath,1,nlev)
c      write(*,*)'compute_columns: nlev=',nlev
      do klev=3,nlev
         if(z(klev).gt.zmin) go to 777
      end do
      klev=nlev
      write(6,*) 'Warning: zmin exceeds z(nlev)',zmin,z(klev)
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
