      subroutine readmodFC(modname,z,t,p,d,h2ovmr,w,roc,nlev,ptrop)
c
c  Reads a model (P,T) from file MODNAME at arbitrary vertical levels and then
c  uses the hydrostatic equation together with the mean molecular weights W to
c  interpolate the T's & P's onto the geometric altitudes specified in Z(NLEV).
c
c  This version adapted to read FASCODE-format model input files
c  Dave Griffith, Sept 00
c
c  Inputs:
c      MODNAME   C*(*)   Name of model
c      Z(NLEV)   R*4     Array of altitudes onto which T & P are to be interpolated
c      W(NLEV)   R*4     Array of mean molecular weights
c      NLEV      I*4     Number of atmospheric levels
c
c  Outputs:
c      T(NLEV)   R*4     Array of interpolated temperatures
c      P(NLEV)   R*4     Array of interpolated pressures
c      D(NLEV)   R*4     Array of interpolated densities
c      h2ovmr(nlev) R*4  H2O vmr from model file
c      ROC       R*4     Radius Of Curvature
c
c  All variables beginning with the letter Z are geometric altitudes
c  All variables beginning with the letter H are geopotential altitudes
c
      implicit none
c      
      character modname*(*),string*100,dummy*20
      integer lunr,nlev,i,k,nlhead,ncol,ninlvl,minlvl,ii,nss
      real*4 gas,radius,ecc2,tlat,gs,gravity,pfact,zero,h2ox,
     & h2ovmr(nlev),h2oold,h2onew,inheight,inmw,
     & pold,zold,hold,told,pnew,hk,x,hnew,tnew,avagadro,log1pxox
      parameter (gas=8.31432,lunr=19,avagadro=6.02217e+23,zero=0.0)
      parameter (minlvl=4000)
      real*4 z(nlev),t(nlev),p(nlev),d(nlev),w(nlev),roc
      real*4 inlvl(minlvl),inpress(minlvl),intemp(minlvl),
     & inh2ovmr(minlvl) !temp storage
      real*8 ptrop
c
c     Read in input levels, pressures and temperatures
	open(unit=lunr,file=modname,status='old')
      read(lunr,*)nlhead,ncol
c      write(*,*)nlhead,ncol
      read(lunr,'(a)')string
c
      if(index(modname,'.mod').gt.0)then           !DG Jan03
c        GFIT format model, listed bottom-up
         call substr(string,dummy,1,nss)  ! nss = number of sub-strings
         if(nss.eq.6) then
            read(string,*)radius,ecc2,tlat,gs,zold,pfact
            ptrop=0.0d0
         elseif(nss.eq.7) then
            read(string,*)radius,ecc2,tlat,gs,zold,pfact,ptrop
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
               read(string,*,end=10)inpress(i),intemp(i),inheight,inmw,
     &         inh2ovmr(i)
            elseif(nss.eq.4) then
               read(string,*,end=10)inpress(i),intemp(i),inheight,inmw
               inh2ovmr(i)=0.0
            elseif(nss.eq.3) then
               read(string,*,end=10)inpress(i),intemp(i),inheight
               inh2ovmr(i)=0.0
            elseif(nss.eq.2) then
               read(string,*,end=10)inpress(i),intemp(i)
               inh2ovmr(i)=0.0
            else
               write(*,*)'READMODFC: Unknown model format: ',modname,nss
               write(*,*)'READMODFC: Last line read: ',string
               stop
            endif
            if(inpress(i).le.0.0)goto 10
         enddo
         stop 'increase parameter minlvl'
10       ninlvl=i-1
c         write(*,*)'ninlvl=',ninlvl
         close(lunr)
      elseif(index(modname,'.zpt').gt.0)then           !DG Jan03
c        FASCOD format model, levels listed top-down
         ninlvl=ncol
         read(lunr,*)(inlvl(i),i=1,ninlvl)
         read(lunr,*)
         read(lunr,*)(inpress(ninlvl+1-i),i=1,ninlvl)
         do i=1,ninlvl
            inpress(i)=inpress(i)/1013.25
         enddo
         read(lunr,*)
         read(lunr,*)(intemp(ninlvl+1-i),i=1,ninlvl)
         radius=6378.00
         ecc2=6.0e-5
         tlat=-34.0
         gs=9.81
         zold=0.03
         pfact=1.0
         close(lunr)
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
      h2oold=inh2ovmr(1)

      pnew=inpress(2)
      tnew=intemp(2)
      h2onew=inh2ovmr(2)
      ii=2
c
      x=tnew/told-1
      hnew=hold+gas*told*log(pold/pnew)/gs/w(1)/log1pxox(x)
c      write(*,*)hnew,hnew/(1-hnew/radius),pnew,tnew
c      if(z(1).lt.zold)
c     &   write(6,*)' Warning! Levels may not extend low enough'
      do k=1,nlev         ! loop over levels
        hk=z(k)/(1+z(k)/radius) ! Convert from geometric to geopotential
 4      if(hk.gt.hnew) then ! hold & hnew are both below hk; read another record
          pold=pnew
          told=tnew
          hold=hnew
          h2oold=h2onew 
c          read(lunr,*,end=3)pnew,tnew
          ii=ii+1
          if(ii.gt.ninlvl)goto 3
          pnew=inpress(ii)
          tnew=intemp(ii)
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
          x=(tnew/told-1)*(hk-hold)/(hnew-hold)
          t(k)=told*(1+x)
          p(k)=pold/pfact*exp((hold-hk)*w(k)*gs*log1pxox(x)/gas/told)
c  density units are cm-2.km-1. multiply density by 10e-5 to get units of
c  cm-3. multiply density by 10e+0 to get units of m-3.
c          d(k)=10132.5*avagadro*p(k)/t(k)/gas  !units of cm-2.km-1
          d(k)=0.101325*avagadro*p(k)/t(k)/gas  !units of cm-3
           h2ox=(h2onew/h2oold-1)*(hk-hold)/(hnew-hold) !RAW Linear
           h2ovmr(k)=h2oold*(1+h2ox) !RAW Linear interpolation of H2O
        endif
c       write(*,*)zold,pold,z(k),t(k),p(k)
      end do
      return
c
3     write(*,*)'modname z(k) zold znew = ',modname,ninlvl,z(k),
     & hold/(1-hold/radius),hnew/(1-hnew/radius)
      stop ' Model levels do not extend high enough'
      end
