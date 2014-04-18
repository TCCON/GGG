      subroutine read_model_fc
     & (lunr,modname,z,w,nlev,
     & t,p,d,h2ovmr,roc,ztrop_ncep,ztrop_gct,tlat)
c
c  Reads a model (P,T) from file MODNAME at arbitrary vertical levels and then
c  uses the hydrostatic equation together with the mean molecular weights W to
c  interpolate the T's & P's onto the geometric altitudes specified in Z(NLEV).
c
c  This version adapted to read FASCODE-format model input files
c  Dave Griffith, Sept 00
c
c  Inputs:
c      LUN_MOD    I*4     Logical Unit number
c      MODNAME    C*(*)   Name of model
c      Z(NLEV)    R*4     Array of altitudes onto which T & P are to be interpolated
c      W(NLEV)    R*4     Array of mean molecular weights
c      NLEV       I*4     Number of atmospheric levels
c
c  Outputs:
c      T(NLEV)    R*4     Array of interpolated temperatures
c      P(NLEV)    R*4     Array of interpolated pressures
c      D(NLEV)    R*4     Array of interpolated densities
c      h2ovmr(nlev) R*4  H2O vmr from model file
c      ROC        R*4     Radius Of Curvature
c      ztrop_ncep R*8    Tropopause altitude from NCEP
c      ztrop_gct  R*8     Tropopause altitude internally computed
c      TLAT       R*4     Latitude
c
c  All variables beginning with the letter Z are geometric altitudes
c  All variables beginning with the letter H are geopotential altitudes
c
      implicit none
      include "../ggg_const_params.f"
      include "const_params.f"
c      
      character modname*(*),string*100,dummy*20
      integer lunr,nlev,i,k,nlhead,ncol,ninlvl,minlvl,ii,nss
      real*4 
     & radius,ecc2,gs,gravity,pfact,
     & h2ox,tlat,
     & h2ovmr(nlev),h2oold,h2onew,inmw,iztrop,
     & pold,zold,hold,told,pnew,hk,x,znew,hnew,tnew,
     & log1pxox
      parameter (minlvl=4000)
      real*4 z(nlev),t(nlev),p(nlev),d(nlev),w(nlev),roc
      real*4 inlvl(minlvl),inpress(minlvl),intemp(minlvl),
     & inh2ovmr(minlvl),inheight(minlvl) !temp storage
      real*8 ztrop_ncep,ptrop_ncep,fr,ztrop_gct,lr,lrwas,zlr,zlrwas
c
      lr=0.0d0 ! avoid compiler warning (may be used uninitialized)
      zlr=0.0d0 ! avoid compiler warning (may be used uninitialized)
c      write(*,'(a,a)')' readmodFC: modname = ',modname
      radius=6378.137  ! Equatorial radius (km)
      ztrop_gct=0.0  ! initial values
      ztrop_ncep=0.0 ! initial values
      ecc2=6.0e-5
      tlat=-34.0
      gs=9.81
      zold=0.0
      pfact=1.0
      ptrop_ncep=200.0d0  ! mbar  Default value in case missing from .mod file
c
      if(index(modname,'mod').gt.0)then           !DG Jan03
c        Read in input levels, pressures and temperatures
         open(unit=lunr,file=modname,status='old')
         read(lunr,*)nlhead,ncol
c         write(*,*)nlhead,ncol
         read(lunr,'(a)')string
c        GFIT format model, listed bottom-up
         call substr(string,dummy,1,nss)  ! nss = number of sub-strings
         if(nss.eq.6) then
            read(string,*)radius,ecc2,tlat,gs,zold,pfact
            ptrop_ncep=0.0d0
         elseif(nss.eq.7) then
            read(string,*)radius,ecc2,tlat,gs,zold,pfact,ptrop_ncep
         else
            write(*,*) 'nss=',nss
            write(*,*) string
            stop 'readmodFC: Unrecognized model header'
         endif

         do k=3,nlhead
            read(lunr,*)
         enddo

         do i=1,minlvl
            read(lunr,'(a)',end=10) string
            call substr(string,dummy,1,nss)  ! nss = number of sub-strings
            if(nss.eq.5) then
             read(string,*,end=10)inpress(i),intemp(i),inheight(i),inmw,
     &         inh2ovmr(i)
            elseif(nss.eq.4) then
              read(string,*,end=10)inpress(i),intemp(i),inheight(i),inmw
               inh2ovmr(i)=0.0
            elseif(nss.eq.3) then
               read(string,*,end=10)inpress(i),intemp(i),inheight(i)
               inh2ovmr(i)=0.0
            elseif(nss.eq.2) then
               read(string,*,end=10)inpress(i),intemp(i)
               inh2ovmr(i)=0.0
               inheight(i)=i
            else
               write(*,*)'READMODFC: Unknown model format: ',modname,nss
               write(*,*)'READMODFC: Last line read: ',string
               stop
            endif
            if(inpress(i).le.0.0)goto 10
         enddo
         stop 'increase parameter minlvl'
10       ninlvl=i-1
c        write(*,*)'minlvl=',ninlvl
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
         zold=inlvl(1)
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
            inh2ovmr(i)=0.0
         enddo
         read(lunr,*)            ! blank
         read(lunr,'(a)')string  ! temperature
         read(lunr,*)            ! $
         read(lunr,'(5(f10.2,1x))')(intemp(i),i=1,ninlvl)
         close(lunr)
         zold=inlvl(1)
      else  
         write(*,*)'Model format not recognised'
         write(*,*)'Model name: ',modname
         stop
      endif
c     Begin
c      if(radius.le.6300.0) write(6,*)'wrong radius:',radius
      if(gs.gt.9.5) gs=gravity(tlat,zero)      ! Surface Gravity (Earth)
      hold=zold/(1+zold/radius)  ! Convert geometric to geopotential altitude
      roc=radius                 ! don't correct for Earth's ecentricity
c                                ! OK for ground based & mid latitude
c      read(lunr,*)pold,told
c      read(lunr,*,end=3)pnew,tnew     
      pold=inpress(1)
      told=intemp(1)
      zold=inheight(1)
      h2oold=inh2ovmr(1)
c      write(*,*) 'nlev=',nlev,pold,told,zold

      if(nlev.gt.1) then
         pnew=inpress(2)
         tnew=intemp(2)
         znew=inheight(2)
         h2onew=inh2ovmr(2)
         lr=-5.0
         zlr=0.0
         ii=2
      else
         hnew=1.E+36
         pnew=pold
         tnew=told
         znew=zold
         h2onew=h2oold
      endif
c     ztrop_gct=10.0
c
c      write(*,*)hold,pold,told
      x=tnew/told-1
      hnew=hold+gas*told*log(pold/pnew)/gs/w(1)/log1pxox(x)
c      write(*,*)'hnew=',hnew,pold,pnew,told,tnew,x
c      write(*,*)'hnew=',hnew,hnew/(1-hnew/radius),pnew,tnew
c      if(z(1).lt.zold)
c     &   write(6,*)' Warning! Levels may not extend low enough'
      iztrop=0.0
      do k=1,nlev                ! loop over levels
        hk=z(k)/(1+z(k)/radius)  ! Convert from geometric to geopotential
 4      if(hk.gt.hnew) then      ! hold & hnew are both below hk; read another record
          pold=pnew
          told=tnew
          zold=znew
          hold=hnew
          lrwas=lr
          zlrwas=zlr
          h2oold=h2onew 
c          read(lunr,*,end=3)pnew,tnew
          ii=ii+1
          if(ii.gt.ninlvl) goto 3
          pnew=inpress(ii)
          tnew=intemp(ii)
          znew=inheight(ii)
c
c  Compute tropopause altitude
          lr=(tnew-told)/(znew-zold) ! Lapse Rate
          zlr=0.5*(zold+znew)        ! Altitude at which lapse rate = lr
c      Find first instance of lapse-rate exceeding -2K/km
          if(abs(radius-6378).lt.50) then     ! Earth
c             write(*,*) k, zold, told, zlr, lr
c            if(ztrop_gct.eq.0.0 .and. zold.gt.5 .and. lr.gt.-2.0) then
             if(iztrop.eq.0.0 .and. zold.gt.5 .and. lr.gt.-2.0) then
                iztrop=1.0 ! do not compute the trop altitude more than once
                ztrop_gct=zlrwas+(zlr-zlrwas)*(-2-lrwas)/(lr-lrwas)
                ztrop_gct=ztrop_gct/(1-ztrop_gct/radius)  ! convert H to Z
             endif
          else                        ! Mars
             ztrop_gct=10.0            ! Mars
          endif
c          
          if(tnew.lt.0.0) then 
             write(*,*) modname
             write(*,*) pnew,tnew
             stop ' Model temperatures must be in Kelvin'
          endif
          x=tnew/told-1
          hnew=hold+gas*told*log(pold/pnew)/gs/w(k)/log1pxox(x)
c      write(*,*)hnew,hnew/(1-hnew/radius),pnew,tnew
          h2onew=inh2ovmr(ii)
          go to 4
        else  !   z(k) is bracketed by zold & znew; proceed with interpolation
           h2ox=0.0
           if(hnew.eq.hold) then
              x=0.0
           else
              x=(tnew/told-1)*(hk-hold)/(hnew-hold)
              if(h2oold.gt.0.0) then
                h2ox=(h2onew/h2oold-1)*(hk-hold)/(hnew-hold)
              endif
           endif
           t(k)=told*(1+x)
           p(k)=pold/pfact*exp((hold-hk)*w(k)*gs*log1pxox(x)/gas/told)
c           write(*,*)'p(k)=', p(k),told,tnew,hk,hold,hnew,gas,w(k),gs,x
           d(k)=0.101325*avagadro*p(k)/t(k)/gas  !units of cm-3
           h2ovmr(k)=h2oold*(1+h2ox) !RAW Linear interpolation of H2O
        endif
c       write(*,*)zold,pold,told,z(k),t(k),p(k),h2ovmr(k)
      end do

c  Convert NCEP tropopause pressure to geometric altitude
      if(ptrop_ncep.gt.0 .and. nlev.gt.1) then
         do k=2,nlev                ! loop over levels
           if(p(k).lt.ptrop_ncep/pfact) exit
         end do
         fr=log(ptrop_ncep/pfact/p(k))/log(p(k-1)/p(k))
         ztrop_ncep=fr*z(k-1)+(1-fr)*z(k)
      else
         ztrop_ncep=0.0  ! 
      endif
c      write(*,'(a40,2f9.2)')modname(:40),ztrop_ncep,ztrop_gct
c      write(51,'(a40,3f9.2)')modname(:40),
c     & ptrop_ncep,ztrop_ncep,ztrop_gct

      return
c
3     write(*,*)'ii, ninlvl, modname z(k) zold znew = ',ii,ninlvl,
     & modname,ninlvl,z(k),
     & hold/(1-hold/radius),hnew/(1-hnew/radius)
      write(*,*)' Warning: Model levels do not extend high enough'
      stop ' Model levels do not extend high enough'
      end
