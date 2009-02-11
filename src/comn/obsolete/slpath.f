      subroutine slpath(nlev,z,t,p,aszain,fovr,roc,zobs,wavtkr,
     & wavmic,zmin,bend,sp,ifail)
c  calculates the geometrical slant paths sp(nlev), averaged over a circular
c  field of view of radius fovr degrees, for a ray of zenith
c  angle aszain degrees, hitting an observer at an altitude hobs km,
c  at the lowest point along the ray trajectory the latitude is tlat degrees
c  and the ray direction has a bearing azim degrees.
c
c  the atmospheric temperature and pressure t(nlev) kelvins & p(nlev) atmospheres
c  is tabulated at altitudes z(nlev) km, which need not be equally spaced.
c
c  input parameters:
c    nlev      integer*4   number of model levels used to discretize atmosphere
c    z(nlev)      real*4   altitudes of model levels
c    t(nlev)      real*4   temperatures of model levels
c    p(nlev)      real*4   model level pressures.
c    aszain(+ve)  real*4   astronomical solar zenith angle (unrefracted).
c    aszain(-ve)  real*4   apparent solar zenith angle (refracted).
c    fovr         real*4   angular radius of instrumental field of view.
c    hobs         real*4   altitude of instrument above sea level.
c    wavtkr       real*4   mean wavenumber of suntracker sensor sensitivity.
c    wavmic       real*4   central wavenumber of microwindow.
c
c  output parameters:
c    hmin         real*4   minimum altitude encountered along ray path.
c    bend         real*4   angle through which ray is bend by refraction.
c    sp(nlev)     real*4   total geometrical slant paths.
c    ifail        integer*4   error flag:
c      ifail=0  sucessful execution
c      ifail=1  apparent solar angle didn't converge within maxiter iterations
c      ifail=2  ray did not intersect atmosphere. sp(j) all = 0.0
c      ifail=3  did not encounter tangent point or z(nlev) within maxstep steps
c      ifail=4  unrefracted ray passed more than 100km below surface.
c      ifail=5  iteration to estimate apparent angle did not converge.
      implicit none
      integer*4 maxiter,ifail,lev1,lev2,nlev,it,iter,k,l,jt,mlev
      parameter (mlev=250)
      real*4 z(nlev),p(nlev),t(nlev),sp(nlev),spg(mlev),aszain,
     & asza,fovr,
     & zobs,wavtkr,wavmic,zmin,zmin0,bend,bend0,bendx,rsza,del_zlev,
     & pobs,tobs,con,ztan,ptan,ttan,prdiff,
     & del_ztan,del_rsza,dum,xsza,weight,rsh,xx,dbdt,b,roc
      real*8 q,opcon,pi,d2r,ri,rizobs,riztan,droc,rtnt
      parameter (pi=3.1415926536d0,d2r=pi/180.d0,maxiter=22)
      ifail=0
      q=opcon(wavtkr)
      rsza=abs(aszain)
      droc=dble(roc)
      if(nlev.gt.mlev) stop 'nlev>mlev'
      if(aszain.lt.0.0) go to 7   ! angle is already refracted - no iteration.
c      write(6,*)'in slpath: wavmic: ',wavmic
c=============================================================================
c  estimate, using a simple empirical model that involves no ray tracing,
c  the refracted angle, rsza, that gives rise to the astronomical angle = asza
c  the purpose of this is to reduce the number of iterations that the ray
c  has to be traced through the atmosphere. this is achieved not just by a
c  more accurate starting guess, but also better partial differentials, dbdt.
c  tests have shown that the rsza estimate is good to better than 0.02 degrees
c  with any model for both ground-based and atmos geometries.
      asza=rsza
      call find(nlev,z,zobs,lev1,lev2)
      del_zlev=z(lev1)-z(lev2)
      rsh=dlog(dble(p(lev2)/p(lev1)))/del_zlev
      pobs=p(lev1)*exp((z(lev1)-zobs)*rsh)
      tobs=t(lev1) - (z(lev1)-zobs)*(t(lev1)-t(lev2))/del_zlev
      rizobs=ri(tobs,pobs,q)
      con=.035
      do 5 it=1,maxiter
      xx=1.0/(con+abs(cos(rsza*d2r)))
      bend=194.0*pobs*sin(rsza*d2r)*con*xx/tobs
      dbdt=1.25*bend*xx
      if(rsza.gt.90.0) then
c         use snell's law in spherical geometry to evaluate true tangent ht.
          rtnt=rizobs*(droc+zobs)*sin(rsza*d2r)
          ztan=rtnt-droc
          if(ztan.lt.z(1)-99.999) then
             ifail=4
             return
          endif
          do 4 jt=1,maxiter
          call find(nlev,z,ztan,lev1,lev2)
          del_zlev=z(lev1)-z(lev2)
          rsh=dlog(dble(p(lev2)/p(lev1)))/del_zlev
          ptan=p(lev1)*exp((z(lev1)-ztan)*rsh)
          ttan = t(lev1) - (z(lev1)-ztan)*(t(lev1)-t(lev2))/del_zlev
          riztan=ri(ttan,ptan,q)
          prdiff=riztan+(droc+ztan)*(1.0-riztan)*rsh
c      write(1,*)ztan,riztan,prdiff,del
          if(abs(prdiff).lt.0.66) go to 8
          del_ztan = (rtnt-(droc+ztan)*riztan)/prdiff
          ztan=ztan+del_ztan
          if(abs(del_ztan).lt.0.001) go to 8
4         continue
c         assumes ray bending is directly proportional to tangent density.
8         b=388.0*ptan/ttan
          bend=-bend+b
          dbdt=-dbdt+1.72*b*(droc+zobs)*rsh/xx
      endif
      del_rsza=asza-rsza-bend
      if(rsza.gt.25.0) del_rsza=rsza*del_rsza/(rsza+dbdt)
      rsza=rsza+del_rsza
      if(abs(del_rsza).lt.0.0003) go to 6
5     continue
      ifail=5
c
c  iterate to find the apparent (refracted) angle, rsza, such that rays
c  of wavenumber wavtkr, when traced back through the atmosphere, emerge
c  at the desired astronomical angle asza. the estimate of the previous
c  program section is used as a starting guess for the iteration. speed of
c  convergence is optimised by using partial differentials dbdt calculated
c  in the previous program section.
6     do 11 iter=1,2*maxiter
      call getsp(nlev,z,p,t,rsza,roc,zobs,zmin0,q,spg,bend,ifail)
      if(rsza.gt.90.0) then
      call getsp(nlev,z,p,t,180.0-rsza,roc,zobs,dum,q,spg,bendx,ifail)
          bend=2.0*bend+bendx
      endif
      del_rsza=asza-rsza-bend
      if(rsza.gt.25.0) del_rsza=rsza*del_rsza/(rsza+dbdt)
c      write(*,*)iter,rsza,bend,zmin0,del_rsza
      rsza=rsza+del_rsza
      if(abs(del_rsza) .lt. 0.00004) go to 12
11    continue
      ifail=1
12    continue
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
      q=opcon(wavmic)
      do l=-1,1
         weight=0.750-0.625*l*l
         xsza=rsza+l*fovr
         zmin=zmin0     ! we want the value when l=0, not l=1
         bend=bend0
         call getsp(nlev,z,p,t,xsza,roc,zobs,zmin0,q,spg,bend0,ifail)
c         write(*,*)'called getsp',xsza,ifail,bend0,sp(1),sp(2),sp(3)

         if(xsza.gt.90.0) then
             do k=1,nlev
                sp(k)=sp(k)+2*weight*spg(k)
             enddo
             call getsp(nlev,z,p,t,180.0-xsza,roc,zobs,dum,q,spg,bendx,
     &       ifail)
             bend0=2*bend0+bendx
         endif
         do k=1,nlev
            sp(k)=sp(k)+weight*spg(k)
         end do
      end do  ! l=-1,1
      return
      end
c
      subroutine getsp(nlev,z,p,t,ths,roc,zobs,zmin,q,spg,bend,
     &ifail)
c  calculates the geometrical slant paths sp(nlev) for a ray of apparent
c  zenith angle ths degrees, hitting an observer at an radius zobs km,
c  from the centre of curvature.
c
c  the atmospheric temperature and pressure t(nlev) & p(nlev) is tabulated
c  at radii z(nlev), which need not be equally spaced.
c
c  the subroutine also returns bend, the angle (degrees) through which the
c  ray was deviated by refraction. this permits getsp to be put into an
c  iteration loop tp find the angle ths which gives rise to a particular
c  astronomical solar zenith angle.
      implicit none
      integer*2 iflag,maxstep
      integer*4 ifail,lev1,lev2,nlev,k
      real*4 z(nlev),p(nlev),t(nlev),spg(nlev),zobs,zmin,bend,ths,
     & ztoa,zhp,del_zlev,pres,temp,x1,x2,a1,a2,con,roc
      real*8 q,ds,dx,dz,dt,dphi,phi,th,thp,zh,rsh,rp,grad,d2r,pi,piby2,
     & droc
      parameter (maxstep=3333)
c
      droc=dble(roc)
      piby2=dacos(0.0d0)
      pi=2*piby2
      d2r=pi/180.0d0
c
      do k=1,nlev
         spg(k)=0.0
      end do
      iflag=0
      zmin=zobs
      bend=0.0
      ztoa=z(nlev)
c
c  if observer is outside atmosphere, start integration at top of atmosphere
      th=dble(ths)*d2r
      if(zobs.gt.ztoa) then
          if(th.lt.piby2) return
          zh=dble(ztoa)
c          zmin=dsin(th)*(droc+zobs)-droc
          zmin=droc*(dsin(th)-1.0d0) + dsin(th)*zobs
          if( zmin.gt.ztoa ) then
              ifail=2  !  ray did not enter atmosphere
              return
          endif
          th=pi-dasin(dsin(th)*(droc+zobs)/(droc+zh))
          phi=ths*d2r-th
      else
          zh=dble(zobs)
          phi=0.0d0
      endif
c
c  ds = integration step  < (layer thickness)/2, or 2km horizontal.
      call find(nlev,z,sngl(zh),lev1,lev2)
      dx=0.5               ! max allowed horizontal step size
      dz=(z(lev1)-z(lev2))/8   ! max allowed vertical step size
      ds=1.0d0/dmax1(1/dx,dabs(dcos(th)/dz))
      dz=ds*dcos(th)
      dx=ds*dsin(th)
c      write(*,*)'getsp:',ds,dz,dx,dcos(th)
      zhp=zh+dz/2.
      dphi=dx/(droc+zhp)
      call find(nlev,z,zhp,lev1,lev2)
      del_zlev=z(lev1)-z(lev2)
      temp=t(lev1) - (z(lev1)-zhp)*(t(lev1)-t(lev2))/del_zlev
      rsh=dlog(dble(p(lev2)/p(lev1)))/del_zlev
      pres=p(lev1)*exp((z(lev1)-zhp)*rsh)
      dt=dx*grad(temp,pres,rsh,q)-dphi
      thp=th+dt/2.
c
c  hp & thp are the altitude and local zenith angle estimated for the middle
c  of the next path segment.
      do k=1,maxstep
         dz=ds*dcos(thp)
         dx=ds*dsin(thp)
         dphi=dx/(droc+zhp)
         call find(nlev,z,zhp,lev1,lev2)
         del_zlev=z(lev1)-z(lev2)
         temp=t(lev1) - (z(lev1)-zhp)*(t(lev1)-t(lev2))/del_zlev
c
c  pressure is expressed as a linear combination of pressures at levels lev1 & lev2.
c  the coefficients a1 & a2 are exponential functions of height, so a1+a2 < 1.
c  a1 & a2 are also the contributions of the current path segment to the slant
c  paths of levels lev1 & lev2 respectively.
         x1=(z(lev1)-zhp)/del_zlev
         x2=x1-1.0
         rp=dlog(dble(p(lev2)/p(lev1)))
         a1=-(x2*dexp(x1*rp))
         a2=x1*dexp(x2*rp)
         pres = p(lev2)*a2 + p(lev1)*a1
         dt=dx*grad(temp,pres,rp/del_zlev,q)-dphi
c
         spg(lev1)=spg(lev1)+a1*ds
         spg(lev2)=spg(lev2)+a2*ds
         phi=phi+dphi
         zh=zh+dz
         th=th+dt
c        write(*,*)zh,th,phi,dz
c
c  terminate integration cleanly if a tangent point is detected.
         if(iflag.eq.1) go to 25
         if( ths.gt.90.0d0 ) then   ! if the initial angle was > 90
            if( (th-piby2)*(th-piby2-dt).le.0.0d0) then
                iflag=1
                con=(piby2-th)/dt
                dz=dz*con
                dt=dt*con
                ds=ds*con
            endif
         endif
c
         zhp=zh+dz/2.
         thp=th+dt/2.
         if(zhp.gt.ztoa) go to 25
      end do ! k=1,maxstep
      ifail=3
c
25    bend=(th+phi)/d2r-dble(ths)
      zmin=min(zobs,sngl(zh))
      return
      end

      subroutine find(nlev,z,zhp,lev1,lev2)
c  finds the model levels which bracket the altitude zhp.
      integer*4 nlev,lev1,lev2
      real*4 z(nlev),zhp
c
      do lev2=1,nlev
         if(z(lev2).ge.zhp) go to 11
      end do
c      write(*,*)'find: ray above top model level',zhp,z(nlev)
      lev2=nlev
c
11    if(lev2.eq.1)then
c       write(*,*)'find: ray below lowest model level',zhp,z(1)
        lev2=2
        endif
      lev1=lev2-1
      return
      end

c
      function grad(temp,pres,rsh,q)
c  calculates  1/n.dn/dz  where n=refractive index.
c  rsh = reciprocal scale height
      real*4 temp,pres
      real*8 q,d,grad,rsh
      d=353.d0*pres*q/temp
c  refractive index = sqrt(1.0d0+2.0d0*d)/(1.0d0-d))
      grad=1.5*d*rsh/((1.0d0+2.0d0*d)*(1.0d0-d))
      return
      end

      function opcon(wav)
      real*4 wav
      real*8 opcon,aa,a,b,c,xt
c  precompute optical constant used to calculate refractive index gradient
      xt=wav*wav*1.0e-8
      a=64.328d0
      b=29498.1d0/(146.0d0-xt)
      c=255.4d0/(41.0d0-xt)
      aa=1.0d-6*(a+b+c)
      opcon=aa*(aa+2.0d0)/(1.225014d0*(aa*(aa+2.0d0)+3.0d0))
c      opcon=0.0  !  disables refraction
      return
      end
c

      function ri(t,p,q)
      real*4 t,p
      real*8 d,q,ri
      d=353.0d0*p*q/t
      if(d.gt.0.97d0)then
         ri=9.9d0
      else
         ri=dsqrt((1.0d0+2.0d0*d)/(1.0d0-d))
      endif
      return
      end
