      subroutine read_model_fc (lunr,modname,z,w,nlev,
     & t,p,d,h2o_dmf,radius,tlat)
c
c  Subroutine has two distinct functions:
c  1) Read the contents of the .mod file into the inxxxx vectors.
c  2) Interpolate vectors onto a user-prescribed altitude grid, Z.
c
c  Reads a model (P,T,Z) from file MODNAME sampled on arbitrary
c  vertical levels and then interpolates the T/P profiles onto
c  the geometric altitudes specified in Z(NLEV).
c
c  If the input model does not contain any altitude information
c  (e.g. it is a sonde profile) then the hydrostatic equation is
c  integrated to assign altitudes to the input levels. For this
c  reason only, the MMW and ROC inputs are needed.
c
c  This version adapted to read FASCODE-format model input files
c  Dave Griffith, Sept 00
c
c  Inputs:
c      LUN_MOD       I*4    Logical Unit Number
c      MODNAME       C*(*)  Name of model
c      Z(NLEV)       R*4    Vector of geometric oputput altitudes
c                           onto which T/P are to be regrided
c      W(NLEV)       R*4    Vector of mean molecular weights
c      NLEV          I*4    Number of atmospheric levels
c
c  Outputs:
c      T(NLEV)       R*4    Regridded temperatures
c      P(NLEV)       R*4    Regridded pressures
c      D(NLEV)       R*4    Densities
c      h2o_dmf(nlev) R*4    Regridded H2O vmr profile
c      RADIUS        R*4    Radius Of Curvature
c      TLAT          R*4    Latitude
c
c  Variables beginning with the letter Z are geometrical altitudes
c  Variables beginning with the letter H are geopotential altitudes
c
c Theoretical Basis:
c  Geopotential altitude (h) is related to geometric altitude (z)
c  by an object of mass m having the same potential energy (PE)
c  when gravity is assumed constant at the standard surface value
c  of 0.980655 m/s/s and therefore independent of altitude and
c  latitude.
c     PE = Integral [m.g(z).dz] = Integral [m.g(0).dh]           (1)
c  where g(z) is the gravitational acceleration at altitude z
c  g(z) = g(0)/[(1+z/R)^2]
c       m.g(0).[R/(1+z/R)-R/(1+0)] = m.g(0).(h-0)
c        h = z/(1+z/R)                                           (2)
c  h is always smaller than z, because if the gravity were
c  constant at the surface value, rather than decreasing with
c  altitude, the object gains PE more quickly with increasing
c  altitude. h and z can be simply inter-converted: equation (2)
c  can be turned inside out
c        z = h/(1-h/R)                                           (3)
c
c  Hydrostatic Equation:  the weight of a thin, unit-area, slab
c  of air matches the pressure difference between the lower and
c  upper surfaces
c      dP = -m(z).g(z).d(z).dz                                   (4)
c      dP = -m(h).g(0).d(h).dh                                   (5)
c  where d(z) is the number density at altitude z
c  where m(z) is the mean molecular mass at altitude z
c  where g(z) is the gravity at altitude z
c  This assumes that the slab is not undergoing acceleration.
c  Outside thunderstorms vertical accelerations do not exceed
c  1mm/s/s, which is 1 part in 10,000 of the gravitational
c  acceleration (9.8 m/s/s).
c
c  (Non-)Ideal Gas Law:   PV=ZnkT  so d=n/V=P/ZkT                (6)
c  where Z is the compressibility factor (non-ideal gas) which for
c  air at 1 atm is 0.9999 @ 300K, 0.9992 @ 250K, and 0.9978 @ 200K.
c  At lower pressures Z tends to 1.0000 at all temperatures. Since
c  temperatures of 200K are only found below 0.2 atm in Earth's
c  atmosphere, ignoring the compressibility of air probably makes
c  a difference of only ~2 parts in 10,000 integrated over the
c  entire atmosphere. But for higher pressures or lower temperature
c  conditions, compressibility would need to be included.
c
c  Combining Hydrostatic & Ideal Gas equations:
c             -m(h).g(0).P
c     dP/dh =  -----------                                       (7)
c                 k.T(h)
c
c             -m(h).g(0)
c     dP/P  =  --------- dh                                      (8)
c                k.T(h)
c  If T and m are constant, then this can be integrated to
c
c                     -(h-h').m(h').g(0) 
c    P(h)/P(h') = Exp[------------------]                        (9)
c                           kT(h')
c
c  Assume T is linear with h:  T(h) = T(h')+LR*(h-h') 
c  where LR is the Lapse Rate (deg/km)
c  Let x = LR*(h-h')/T(h') be the fractional change in T
c  between h' and h.  So dx=dh.LR/T(h')
c     T(x) = T(h')*(1+x)
c
c                -m(h).g(0)   
c       dP/P  =  ----------  dx                                 (10)
c                LR.k.(1+x)
c
c  Assume  m(h) = m(h')+dm/dh.(h-h') = m(h')(1+x.dm/dh.T(h')/LR/m(h'))
c  Let b=dm/dh.T(h')/LR/m(h') so that m(h) = m(h').(1+b.x)
c  b.x = (h-h').dm/dh/m(h') is the fractional change in m in 
c  ascending from h'to h.
c  How large is b.x? At the top of the PBL the H2O mole fraction
c  can change from 0 to 2% in 1 km, which means that m(h) changes
c  from 28.964 to 28.964*(1.00-0.02)+18.015*0.02 = 29.627.
c  So dm/m = 0.02*(1-18.015/29.864) = 0.008 km-1
c  At the top of the PBL the levels are only separated by 0.5 km
c  so the maximum possible value of x.b is 0.5*0.008 = 0.004 
c  So the changing MMW of air is a small effect, even in summer,
c  but needs to be included if high accuracy is required.
c
c                -g(0).m(h').[1+b.x]   
c       dP/P  =   ------------------  dx                        (11)
c                     LR.k.(1+x)
c
c  The integral of (1+b.x)/(1+x) is fortunately of closed-form:
c     Integral[(1+b.x)/(1+x)]dx = [(1-b).loge(1+x)+b.x]
c So
c     log(P(h)/P(h'))=-g(0).m(h')/k/LR.[(1-b).loge(1+x)+b.x]
c  Note that if m(h) is constant, b=0, in which case
c     log(P(h)/P(h'))=-g(0).m(h')/k/LR.loge(1+x]
c                    =-g(0).m(h')/k/LR.loge(1+LR*(h-h')/T(h')]  (12)
c     P(h)/P(h') = [1+LR*(h-h')/T(h')]^[-g(0).m(h)/k/LR]        (13)
c  But if LR=0 a numerical problem arises since P(h)/P(h')=1^Inf
c  So    
c     log(P(h)/P(h')) =-g(0).m(h').(h-h')/k/T(h').[(1-b).loge(1+x)/x+b]
c  which can be rewritten as
c     h-h' = log[P(h)/P(h')].kT(h')/m/g(0)/[b+(1-b).loge(1+x)/x]
c  So
c                      -(h-h').g(0).m(h').[b+(1-b).loge(1+x)/x]
c     P(h)/P(h') = Exp[---------------------------------------] (14)
c                                      kT(h')
c
c  If LR --> 0, x --> 0, loge[1+x]/x --> 1, and so equation (14)
c  reduces to (9). So the use of a robust intrinsic
c  function for loge(1+x)/x avoids the LR=0 problem.
c  Replace the mean molecular mass (m) by the mean molar mass (M)
c  M = A.m, where A is Avocadro's constant. Since The gas constant
c  R=A.k we can write 

c               -R.T(h').loge[P(h)/P(h')]
c     h-h' =  --------------------------------
c             M(h').g(0).[b+(1-b).loge(1+x)/x] 

c  Since b=dM/dh.T(h')/LR/M(h') it can become infinite when LR=0.
c  But in this case x=0 so loge(1+x)/x=1, so the b-terms cancel.
c  So no mathematical problem but still a numberical problem.
c  For small x, loge(1+x)/x ~= 1-x/2
c  and so [b+(1-b).loge(1+x)/x] = [b+(1-b).(1-x/2)] = 1-x/2+b.x/2
c  which is well-behaved because both  b.x = (h-h').(dm/dh)/m(h')
c  and x = (h-h').LR/T(h')  are << 1.
c
c  A bettewr approximation to loge(1+x)/x is the Pade (6+x)/(6+4x)
c  So [b+(1-b).loge(1+x)/x] --> b + (1-b).(6+x)/(6+4x)
c                            =  [b(6+4x)+(1-b).(6+x)]/(6+4x)
c                            =  [6+3b.x+x]/(6+4x)
c  which is numerically stable because the b is multiplied by x.
c
c  An even better Pade approximation to loge(1+x)/x is
c
c               (15-4(x/(2+x))**2)
c    f(x) = ---------------------------
c           (15-9(x/(2+x))**2).(1+x/2)
c
c  This is phenomenally accurate for x > -0.9
c So we can write

c  We want to compute  g(x) = b + (1-b).f(x)
c  Let z=(x/(2+x))**2
c
c          1 +b.x/2 +5z(1-b)/(15-9z)
c  g(x) =  ------------------------
c                  (1+x/2)
c
c  For small x values, x^2 becomes negligible as does z and so
c  g(x) = (1+b.x/2)/(1+x/2) = 1 -x/2 +b.x/2


      implicit none

      real*4 zero
      parameter(zero=0.0)

      include "const_params.f"
c      
      character modname*(*),string*100,dummy*20
      integer lunr,nlev,i,k,nlhead,ncol,ninlvl,minlvl,ii,nss,iflag
      real*4 
     & radius,ecc2,gs,gravity,pfact,
     & tlat,
     & beta,lr,tbar,
     & h2o_dmf(nlev),wvold,wvnew,
     & pold,hold,told,pnew,hk,x,hnew,tnew,
     & log1pxox
      parameter (minlvl=4000)
      real*4 z(nlev),t(nlev),p(nlev),d(nlev),w(nlev)
      real*4 inlvl(minlvl),inpress(minlvl),intemp(minlvl),
     & inh2o_dmf(minlvl),inheight(minlvl),inmmw(minlvl)  !temp storage
      real*8 rh

      radius=6378.137  ! Equatorial radius (in case missing from .mod)
      ecc2=6.0e-5
      tlat=-34.0
      gs=9.80665
      pfact=1.0
c
      if(index(modname,'mod').gt.0) then           !DG Jan03
c        Read in input levels, pressures and temperatures
c         write(*,*) 'modname = ',modname(:70)
         open(unit=lunr,file=modname,status='old')
         read(lunr,*)nlhead,ncol
         read(lunr,'(a)') string  ! radius,ecc2,tlat,gs,hold,pfact
         call substr(string,dummy,1,nss)  ! nss = number of sub-strings
c         write(*,*)' read_model_fc: nss,string = ',nss,string
         if(nss.ge.6) then
            read(string,*)radius,ecc2,tlat,gs,hold,pfact
         else
            write(*,*) 'nss=',nss
            write(*,*) string
            stop 'Line 75: read_model_fc: Unrecognized model header'
         endif

c    Skip column labels and units
         do i=3,nlhead
            read(lunr,'(a)') 
         end do

         do i=1,minlvl
            read(lunr,'(a)',end=10) string
            if(ncol.ge.10) then
               read(string,*,end=10)inpress(i),intemp(i),inheight(i),
     &         inmmw(i),inh2o_dmf(i),rh
            elseif(ncol.eq.6) then
               read(string,*,end=10)inpress(i),intemp(i),inheight(i),
     &         inmmw(i),inh2o_dmf(i),rh
            elseif(ncol.eq.5) then
               read(string,*,end=10)inpress(i),intemp(i),inheight(i),
     &         inmmw(i),inh2o_dmf(i)
            elseif(ncol.eq.4) then
               read(string,*,end=10)inpress(i),intemp(i),inheight(i),
     &         inmmw(i)
               inh2o_dmf(i)=0.0
            elseif(ncol.eq.3) then
               read(string,*,end=10)inpress(i),intemp(i),inheight(i)
               inh2o_dmf(i)=0.0
            elseif(ncol.eq.2) then ! Unreliable or no height info. Use hydrostatic eqtn
               read(string,*,end=10)inpress(i),intemp(i)
               inh2o_dmf(i)=0.0
               inmmw(i)=28.9640
c               inmmw(i)=28.9644
               if(i.le.1) then
                  inheight(i)=hold
               else
                  if(abs(intemp(i-1)-intemp(i)).gt.0.01*intemp(i)) then
                     tbar=(intemp(i-1)-intemp(i))/
     &               log(intemp(i-1)/intemp(i))
                  else
                     tbar=(intemp(i-1)+intemp(i))/2
                  endif
                  inheight(i)=inheight(i-1)+gas*tbar/inmmw(i)/gs*
     &            log(inpress(i-1)/inpress(i))
               endif
c               write(*,*)'i,inheight(i)=',i,inheight(i)
            else
               write(*,*)'read_model_fc: Bad model format:',modname,nss,
     &         ncol
               write(*,*)'read_model_fc: Last line read:',string
               stop
            endif
            if(inpress(i).le.0.0)goto 10
            if(intemp(i).lt.0.0) then 
               write(*,*) modname,i,intemp(i)
               stop ' Model temperatures must be in Kelvin'
            endif
         enddo  !   do i=1,minlvl
         stop 'increase parameter minlvl'
10       ninlvl=i-1
         close(lunr)
      elseif(index(modname,'.zpt').gt.0)then           !DG Jan03
c        FASCOD format model, levels listed top-down
         open(unit=lunr,file=modname,status='old')
         read(lunr,*)nlhead,ncol
         read(lunr,'(a)')string
         ninlvl=ncol
         read(lunr,*)(inlvl(i),i=1,ninlvl)
         read(lunr,*)
         read(lunr,*)(inpress(ninlvl+1-i),i=1,ninlvl)
         do i=1,ninlvl
            inpress(i)=inpress(i)/1013.25
         enddo
         read(lunr,*)
         read(lunr,*)(intemp(ninlvl+1-i),i=1,ninlvl)
         hold=inlvl(1)
         close(lunr)
      elseif(index(modname,'.prf').gt.0)then           ! DG Apr 08
c        FASCOD format model, levels listed top-down
         open(unit=lunr,file=modname,status='old')
         read(lunr,*)            ! blank
         read(lunr,'(a)')string  ! number of levels
         read(lunr,*)            ! $
         read(lunr,*)ninlvl
         read(lunr,*)
         read(lunr,*)
         read(lunr,'(a)')string  ! altitude
         read(lunr,*)            ! $
         read(lunr,'(5(f10.2,1x))')(inlvl(i),i=1,ninlvl)
         read(lunr,*)            ! blank
         read(lunr,'(a)')string  ! pressure
         read(lunr,*)            ! $
         read(lunr,'(5(e10.3,1x))')(inpress(i),i=1,ninlvl)
         do i=1,ninlvl
            inpress(i)=inpress(i)/1013.25
            inh2o_dmf(i)=0.0
         enddo
         read(lunr,*)            ! blank
         read(lunr,'(a)')string  ! temperature
         read(lunr,*)            ! $
         read(lunr,'(5(f10.2,1x))')(intemp(i),i=1,ninlvl)
         close(lunr)
         hold=inlvl(1)
      else  
         write(*,*)'Line 174: Model format not recognised'
         write(*,*)'Model name: ',modname
         stop
      endif


c      write(*,*)'ninlvl=',ninlvl
c  Begin computing T/P at altitudes Z(i)
c      if(radius.le.6300.0) write(6,*)'wrong radius:',radius
      if(gs.gt.9.5) gs=gravity(tlat,zero)      ! Surface Gravity (Earth)
      pold=inpress(1)
      told=intemp(1)
      hold=inheight(1)
      wvold=inh2o_dmf(1)
c      write(*,*) 'nlev=',nlev,pold,told,hold

      if(nlev.gt.1) then
         pnew=inpress(2)
         tnew=intemp(2)
         hnew=inheight(2)
         wvnew=inh2o_dmf(2)
         ii=2
      else
         pnew=pold
         tnew=told
         hnew=hold
         wvnew=wvold
      endif
c
c      write(*,*)hold,pold,told
c  Ignore heights in the .mod file, except for first one.
c  Recompute from pressures and temperatures.
      x=tnew/told-1
      hnew=hold+gas*told*log(pold/pnew)/gs/w(1)/log1pxox(x)
c      write(*,*)' k     hold     told     pold     z(k)     t(k)       
c     & p(k)     lr     beta'
      iflag=0
      ii=1
      do k=1,nlev          ! loop over user-prescribed levels (not model levels)
c ! Convert specified altitude grid from geometric to geopotential
         hk=z(k)/(1+z(k)/radius) 
         do while (inheight(ii+1).lt.hk)
            if(ii.lt.ninlvl) then
               ii=ii+1
            else
               write(*,*)'ii,ninlvl,modname,hk,height = ',ii,ninlvl,
     &         modname,hk,z(k),inheight(ii),inheight(ii+1)
               stop ' Model levels do not extend high enough'
            endif
         end do 
c         write(*,'(a,3(i3,f10.2),f8.3)')
c     & 'k,hk,inheight(ii),ii+1,inheight(ii+1)=',
c     & k,hk,ii,inheight(ii),ii+1,inheight(ii+1)
c    hk is now bracketed by h(ii+1) and h(ii) -- proceed with interpolation
         lr=(intemp(ii+1)-intemp(ii))/(inheight(ii+1)-inheight(ii)) ! lapse rate
         if(lr*(intemp(ii+1)-intemp(ii)).gt.0.01*intemp(ii)) then
            beta=log(1+lr*(hk-inheight(ii))/intemp(ii))/
     &           log(1+lr*(inheight(ii+1)-inheight(ii))/intemp(ii))
         else
            beta=(hk-inheight(ii))/(inheight(ii+1)-inheight(ii))
         endif
         h2o_dmf(k)=inh2o_dmf(ii)+beta*(inh2o_dmf(ii+1)-inh2o_dmf(ii))
         t(k)=intemp(ii)+lr*(hk-inheight(ii))
c         t(k)=intemp(ii)+beta*(intemp(ii+1)-intemp(ii))
         p(k)=inpress(ii)*(inpress(ii+1)/inpress(ii))**beta
         p(k)=p(k)/pfact  ! convert to atmospheres
         d(k)=0.101325*avogadro*p(k)/t(k)/gas  !units of cm-3
c         write(*,'(i3,8f10.4)')k,hold,told,pold,z(k),t(k),p(k),lr,beta
      end do    !  do k=1,nlev        ! loop over output levels

      return
      end
