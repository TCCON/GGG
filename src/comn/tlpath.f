      subroutine tlpath(nlev,z,t,p,asza_in,fovr,roc,zobs,wavtkr,
     & wavmic,zmin,bend,sp,ifail)
c  Calculates the geometrical slant paths sp(nlev), averaged over a circular
c  field of view of radius fovr degrees, for a ray of zenith
c  angle asza_in degrees, hitting an observer at an altitude Zobs km,
c  at the lowest point along the ray trajectory the latitude is tlat degrees
c  and the ray direction has a bearing azim degrees.
c
c  the atmospheric temperature and pressure t(nlev) kelvins & p(nlev) atmospheres
c  is tabulated at altitudes z(nlev) km, which need not be equally spaced.
c
c  Input Parameters:
c    nlev      integer*4   number of model levels used to discretize atmosphere
c    z(nlev)      real*4   altitudes of model levels
c    t(nlev)      real*4   temperatures of model levels
c    p(nlev)      real*4   model level pressures.
c    asza_in(+ve) real*4   astronomical solar zenith angle (unrefracted).
c    asza_in(-ve) real*4   apparent solar zenith angle (refracted).
c    fovr         real*4   angular radius of instrumental field of view.
c    zobs         real*4   geometric altitude of instrument above sea level.
c    wavtkr       real*4   mean wavenumber of suntracker sensor sensitivity.
c    wavmic       real*4   central wavenumber of microwindow.
c
c  Output Parameters:
c    zmin         real*4   minimum geometric altitude encountered along ray path.
c    bend         real*4   angle through which ray is bend by refraction.
c    sp(nlev)     real*4   total geometrical slant paths.
c    ifail        integer*4   error flag:
c      ifail=0  sucessful execution
c      ifail=1  apparent solar angle didn't converge within maxiter iterations
c      ifail=2  ray did not intersect atmosphere. sp(j) all = 0.0
c      ifail=3  did not encounter tangent point or z(nlev) within maxstep steps
c      ifail=4  unrefracted ray passed more than 100km below surface.
c      ifail=5  iteration to estimate apparent angle did not converge.
c      ifail=6  tlpath called with NLEV=1

      implicit none

      integer*4 maxiter,ifail,lev1,lev2,nlev,it,iter,k,l,jt

      real*4 z(nlev),p(nlev),t(nlev),sp(nlev),spg(nlev),asza_in,
     & asza,fovr,sd2r,opcon_tkr,opcon_mic,calc_opcon,
     & zobs,wavtkr,wavmic,zmin,zmin0,bend,bend0,bendx,rsza,del_zlev,
     & pobs,tobs,con,ztan,ptan,ttan,prdiff,
     & del_ztan,del_rsza,dum,xsza,weight,rpsh,xx,dbdt,b,roc

      real*8 dpi,d2r,ri,rizobs,riztan,
     & droc,rtnt

      parameter(dpi=3.14159265359d0)
      parameter (d2r=dpi/180.d0,maxiter=22)

      sd2r=acos(0.0)/90.0
      ifail=0
      dbdt=0.0  ! avoid compiler warning (may be used uninitialized)
      if(nlev.eq.1) then
        write(*,*) 'TLPATH exiting due to NLEV<=1'
        sp(1)=0.0
        bend=0.0
        zmin=z(1)
        ifail=6
        return
      endif
      opcon_tkr=calc_opcon(wavtkr)
      rsza=abs(asza_in)
      droc=dble(roc)
      if(asza_in.lt.0.0) go to 7   ! angle is already refracted - no iteration.
c      write(6,*)'in tlpath: wavmic: ',wavmic
c=============================================================================
c  estimate, using a simple empirical model that involves no ray tracing,
c  the refracted angle, rsza, that gives rise to the astronomical angle = asza
c  the purpose of this is to reduce the number of iterations that the ray
c  has to be traced through the atmosphere. this is achieved not just by a
c  more accurate starting guess, but also better partial differentials, dbdt.
c  tests have shown that the rsza estimate is good to better than 0.02 degrees
c  with any model for both ground-based and occultation geometries.
      asza=rsza
      call find(nlev,z,zobs,lev1,lev2)
      del_zlev=z(lev1)-z(lev2)
      rpsh=log(p(lev2)/p(lev1))/del_zlev
      pobs=p(lev1)*exp((z(lev1)-zobs)*rpsh)
      tobs=t(lev1)-(z(lev1)-zobs)*(t(lev1)-t(lev2))/del_zlev
      rizobs=ri(tobs,pobs,opcon_tkr)
      con=.035
      do it=1,maxiter
         xx=1.0/(con+abs(cos(rsza*sd2r)))
         bend=194.0*pobs*sin(rsza*sd2r)*con*xx/tobs
         dbdt=1.25*bend*xx
         if(rsza.gt.90.0) then
c   use snell's law in spherical geometry to evaluate true tangent ht.
             rtnt=rizobs*(droc+zobs)*dsin(rsza*d2r)
             ztan=sngl(rtnt-droc)
             if(ztan.lt.z(1)-99.999) then
                ifail=4
                return
             endif
             do jt=1,maxiter
                call find(nlev,z,ztan,lev1,lev2)
                del_zlev=z(lev1)-z(lev2)
                rpsh=log(p(lev2)/p(lev1))/del_zlev
                ptan=p(lev1)*exp((z(lev1)-ztan)*rpsh)
                ttan=t(lev1)-(z(lev1)-ztan)*(t(lev1)-t(lev2))/del_zlev
                riztan=ri(ttan,ptan,opcon_tkr)
                prdiff=sngl(riztan+(droc+ztan)*(1.0-riztan)*rpsh)
                if(abs(prdiff).lt.0.66) exit
                del_ztan = sngl(rtnt-(droc+ztan)*riztan)/prdiff
                ztan=ztan+del_ztan
                if(abs(del_ztan).lt.0.001) exit
             end do  !   do jt=1,maxiter
c   assume ray bending is directly proportional to tangent density.
             b=388.0*ptan/ttan
             bend=-bend+b
             dbdt=-dbdt+1.72*b*sngl(droc+zobs)*rpsh/xx
         endif
         del_rsza=asza-rsza-bend
         if(rsza.gt.25.0) del_rsza=rsza*del_rsza/(rsza+dbdt)
c         write(*,*)iter,rsza,bend,zmin0,del_rsza
         rsza=rsza+del_rsza
         if(abs(del_rsza).lt.0.0003) go to 6
      end do   !   do it=1,maxiter
      ifail=5
c
c  iterate to find the apparent (refracted) angle, rsza, such that rays
c  of wavenumber wavtkr, when traced back through the atmosphere, emerge
c  at the desired astronomical angle asza. the estimate of the previous
c  program section is used as a starting guess for the iteration. speed of
c  convergence is optimised by using partial differentials dbdt calculated
c  in the previous program section.
6     do  iter=1,2*maxiter
         call getsp(nlev,z,p,t,rsza,roc,zobs,zmin0,opcon_tkr,spg,bend,
     &    ifail)
         if(rsza.gt.90.0) then
            call getsp(nlev,z,p,t,180.0-rsza,roc,zobs,dum,opcon_tkr,spg,
     &      bendx,ifail)
             bend=2.0*bend+bendx
         endif
         del_rsza=asza-rsza-bend
         if(rsza.gt.25.0) del_rsza=rsza*del_rsza/(rsza+dbdt)
         rsza=rsza+del_rsza
         if(abs(del_rsza) .lt. 0.00004) go to 7
      end do   !   do  iter=1,2*maxiter
      ifail=1
c
7     do k=1,nlev
         sp(k)=0.0
      end do
c
c  having found rsza we now trace rays of wavenumber wavmic from the top, middle
c  and bottom of the field of view and weight their slant contributions by
c  0.125, 0.750 & 0.125 respectively. this corrects exactly for quadratic
c  curvature of the atmospheric airmass as a function of angle over the
c  circular field of view.
      opcon_mic=calc_opcon(wavmic)
      do l=-1,1
         weight=0.750-0.625*l*l
         xsza=rsza+l*fovr
         zmin=zmin0     ! we want the value when l=0, not l=1
         bend=bend0
         call getsp(nlev,z,p,t,xsza,roc,zobs,zmin0,opcon_mic,spg,bend0,
     &   ifail)
c         write(*,*)'called getsp',xsza,ifail,bend0,sp(1),sp(2),sp(3)

         if(xsza.gt.90.0) then
             do k=1,nlev
                sp(k)=sp(k)+2*weight*spg(k)
             enddo
             call getsp(nlev,z,p,t,180.0-xsza,roc,zobs,dum,opcon_mic,
     &       spg,bendx,ifail)
             bend0=2*bend0+bendx
         endif
         do k=1,nlev
            sp(k)=sp(k)+weight*spg(k)
         end do
      end do  ! l=-1,1
      return
      end
c
      subroutine getsp(nlev,z,p,t,ths,roc,zobs,zmin,opcon,spg,bend,
     &ifail)
c  Calculates the geometrical slant paths sp(nlev) for a ray of apparent
c  zenith angle ths degrees, hitting an observer at an radius zobs km,
c  from the centre of curvature.
c
c  The atmospheric temperature and pressure t(nlev) & p(nlev) is tabulated
c  at radii z(nlev), which need not be equally spaced.
c
c  The subroutine also returns bend, the angle (degrees) through which the
c  ray was deviated by refraction. this permits getsp to be put into an
c  iteration loop tp find the angle ths which gives rise to a particular
c  astronomical solar zenith angle.
      implicit none
      real*8 dpi          !PI in double precision
      parameter(dpi=3.14159265359d0)

      integer*4 ifail,lev1,lev2,nlev,k,maxstep
      real*4 z(nlev),p(nlev),t(nlev),spg(nlev),zobs,zmin,bend,ths,
     & zmp,del_zlev,pres,temp,x1,x2,a1,a2,con,roc,rpsh,opcon,grad
      real*8 ds,dx,dz,dt,dphi,phi,th,thmp,zh,rp,d2r,
     & piby2,droc
      parameter (maxstep=6400)

      droc=dble(roc)
      piby2=dpi/2.d0
      d2r=dpi/180.0d0
c
      do k=1,nlev
         spg(k)=0.0
      end do
      zmin=zobs
      bend=0.0
c
c  If observer is outside atmosphere, start integration at top of atmosphere
c  assuming no refraction above Z(nlev)
      th=dble(ths)*d2r
      if(zobs.gt.z(nlev)) then
          if(th.lt.piby2) return
          zh=dble(z(nlev))
          zmin=sngl( dsin(th)*(droc+zobs) - droc )
c          write(*,*)ths,th,dsin(th),zobs,z(nlev),zmin
          if( zmin.gt.z(nlev) ) then
              ifail=2  !  ray did not enter atmosphere
              return
          endif
          th=dpi-dasin(dsin(th)*(droc+zobs)/(droc+zh))
          phi=ths*d2r-th
      else
          zh=dble(zobs)
          phi=0.0d0
      endif
c
c  ds = integration step = (layer thickness)/4, or 0.5 km,
c  whichever is the smaller.  For a 1 km layer spacing,
c  these terms become equal when cos(theta)=0.5, theta=60 deg
      call find(nlev,z,sngl(zh),lev1,lev2)
      ds=1.0d0/dmax1(1/0.5d0,4*dabs(dcos(th)/(z(lev1)-z(lev2))))
      dz=ds*dcos(th)
      dx=ds*dsin(th)
c      write(*,*)'getsp:',ths,ds,dz,dx,dcos(th)
      zmp=sngl(zh+dz/2.)
      dphi=dx/(droc+zmp)
      call find(nlev,z,zmp,lev1,lev2)
      del_zlev=z(lev1)-z(lev2)
      temp=t(lev1)-(z(lev1)-zmp)*(t(lev1)-t(lev2))/del_zlev
      rpsh=log(p(lev2)/p(lev1))/del_zlev
      pres=p(lev1)*exp((z(lev1)-zmp)*rpsh)
      dt=dx*grad(temp,pres,rpsh,opcon)-dphi
      thmp=th+dt/2.
c
c  zmp & thmp are the altitude and local zenith angle estimated for the mid-point
c  of the next path segment.
      do k=1,maxstep
         dz=ds*dcos(thmp)
         dx=ds*dsin(thmp)
         dphi=dx/(droc+zmp)
         del_zlev=z(lev1)-z(lev2)
         temp=t(lev1)-(z(lev1)-zmp)*(t(lev1)-t(lev2))/del_zlev
c
c  Pressure is expressed as a linear combination of pressures at levels lev1 & lev2.
c  the coefficients a1 & a2 are exponential functions of height, so a1+a2 > 1.
c  When zmp=z(lev1), x1=0, x2=-1, a1=1, a2=0, 
c  When zmp=z(lev2), x1=1, x2=0, a1=0, a2=1, 
c  When zmp=[z(lev1)+z(lev2)]/2, x1=0.5, x2=-0.5, a1=0.5*exp(0.5*rp), a2=0.5*exp(-0.5*rp), 
c  a1 & a2 are also the contributions of the current path segment to the slant
c  paths of levels lev1 & lev2 respectively.
         x1=(z(lev1)-zmp)/del_zlev
         x2=x1-1.0
         rp=dlog(dble(p(lev2)/p(lev1)))
         a1=-x2*exp(x1*sngl(rp))  !  -(x1-1).exp(x1*rp)
         a2= x1*exp(x2*sngl(rp))  !       x1.exp(x1*rp).exp(-rp)
         pres = p(lev2)*a2 + p(lev1)*a1

c  No density weighting
         a1=1.-x1
         a2=1.-a1
         pres = p(lev1)*exp(x1*sngl(rp))

         dt=dx*grad(temp,pres,sngl(rp)/del_zlev,opcon)-dphi
c
         spg(lev1)=spg(lev1)+a1*sngl(ds)
         spg(lev2)=spg(lev2)+a2*sngl(ds)
         phi=phi+dphi
         zh=zh+dz
         th=th+dt
c
c  Terminate integration cleanly if a tangent point is detected.
         if( ths.gt.90.0d0 ) then   ! if the initial angle was > 90
            if( (th-piby2)*(th-piby2-dt).le.0.0d0) then  ! Gone past TP
c         Back-track to tangent point
               con=sngl((th-piby2)/dt)
               spg(lev1)=spg(lev1)-a1*con*sngl(ds)
               spg(lev2)=spg(lev2)-a2*con*sngl(ds)
               zh=zh-dz*con
               th=th-dt*con
               phi=phi-dphi*con
               go to 25
            endif
         endif
c
c  Reset ray-tracing when level is crossed.
         if(zh.gt.z(lev2)) then    !   Back-track to z(lev2)
            con=sngl((zh-z(lev2))/dz)
            spg(lev1)=spg(lev1)-a1*con*sngl(ds)
            spg(lev2)=spg(lev2)-a2*con*sngl(ds)
            zh=z(lev2)
            th=th-dt*con
            phi=phi-dphi*con
            if(zh.ge.z(nlev)) go to 25
         endif 

         call find(nlev,z,zmp,lev1,lev2)
         zmp=sngl(zh+dz/2.)
         thmp=th+dt/2.

      end do ! k=1,maxstep
      ifail=3
c
25    bend=sngl((th+phi)/d2r)-ths
      zmin=min(zobs,sngl(zh))
      return
      end

      subroutine find(nlev,z,zz,lev1,lev2)
c  Finds the model levels which bracket the altitude zz.
      integer*4 nlev,lev1,lev2
      real*4 z(nlev),zz
c
      do lev2=1,nlev
         if(z(lev2).ge.zz) go to 11
      end do
c      write(*,*)'find: ray above top model level',zz,z(nlev)
      lev2=nlev
c
11    if(lev2.eq.1) then
c       write(*,*)'find: ray below lowest model level',zz,z(1)
        lev2=2
        endif
      lev1=lev2-1
      return
      end

      function calc_opcon(wav)
      real*4 wav
      real*4 calc_opcon,aa,a,b,c,xt
c  Precompute optical constant used to calculate refractive index gradient
      if(wav.le.0.0) then
         calc_opcon=0.0  !  disables refraction
      else
         xt=wav*wav*1.0e-8
         a=64.328
         b=29498.1/(146.0-xt)
         c=255.4/(41.0-xt)
         aa=1.0e-6*(a+b+c)
         calc_opcon=aa*(aa+2.0)/(1.225014*(aa*(aa+2.0)+3.0))
      endif
      return
      end
c
      function ri(temp,pres,opcon)
      real*4 temp,pres,opcon,d
      real*8 ri
      d=353.0*pres*opcon/temp
      if(d.gt.0.97)then
         ri=9.9d0
      else
         ri=dsqrt((1.0d0+2.0d0*d)/(1.0d0-d))
      endif
      return
      end
c
      function grad(temp,pres,rpsh,opcon)
c  Calculates  1/n.dn/dz  where n=refractive index.
c  rpsh = reciprocal scale height
      real*4 temp,pres,rpsh,opcon,d,grad
      d=353*pres*opcon/temp
c  refractive index = sqrt(1.0+2.0*d)/(1.0-d))
      grad=1.5*d*rpsh/((1.0+2.0*d)*(1.0-d))
      return
      end

