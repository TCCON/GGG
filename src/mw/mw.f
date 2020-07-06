c  Program to assist in the identification of spectral windows.
c
c  Calculates approximate transmittance spectra for a user-defined
c  target gas and for all ther other gases by representing the
c  inhomogeneous atmospheric ray path by a homogenous path defined by
c  its absorber-weighted pressure and temperature, which are calculated
c  separately for each gas.
c
c  For every spectral point, it calculates the
c      target = transmittance due to a particular target gas
c      tother = transmittance due to all other gases
c  It then determines the very best window: highest (1-target)/(1-other+noise)
c  Prints out all the spectral points for which the snr > smax/8.
c
c  Program will require 24 Mbyte of memory per 1000 cm-1 searched.
c
      implicit none
      include "../gfit/ggg_int_params.f"

      integer*4 j,lunr,nlhead,ncol,igas,ngas,ilev,nlev,
     & kgas,kiso,reclen,jline,kline1,kline2,nlines,posnall,molno,iso,
     & lmax,mmp,nmp,ip,lunw,lunt,linelist,kk,i,nfft,idum
      integer*8 fsib,file_size_in_bytes
      parameter (lunr=14,lunw=16,lunt=17)
      parameter (mmp=512*1024)

      real*4 z(mlev),sp(mlev),t(mlev),p(mlev),d(mlev),
     & vmr(mgas,mlev),tfac(mgas),slp,slpv,roc,zobs,zmin,
     & slcol(mgas),target(mmp),tother(mmp),spts(mmp),
     & gse(mmp),vac(mmp,mgas),
     & pbar(mgas),tbar(mgas),wbar(mgas),sxtot(mgas),sig_other,
     & smax,snr,noise,theta,xx,zz,curv_target,curv_other

      real*8 fmin,fmax,dnu,freq,sx,stren,abhw,sbhw,eprime,tdabhw,
     & pshift,d2r,findex,frac,dpi
      parameter(dpi=3.14159265359d0)
      parameter (noise=0.005,dnu=0.01d0,roc=6356.0,d2r=dpi/180.0d0)
c
      character string*1000,llformat*53,rotate*24,llname*80,version*50
      character menuinputfile*40
c
      version=' MW         Version 4.11      2018-12-01      GCT '

      idum=mauxcol     ! Avoid compiler warning: Unused parameter
      idum=mcolvav     ! Avoid compiler warning: Unused parameter
      idum=mfilepath   ! Avoid compiler warning: Unused parameter
      idum=mrow_qc     ! Avoid compiler warning: Unused parameter
      idum=mspeci      ! Avoid compiler warning: Unused parameter
      idum=mvmode      ! Avoid compiler warning: Unused parameter
      idum=ncell       ! Avoid compiler warning: Unused parameter
      idum=nchar       ! Avoid compiler warning: Unused parameter

c  Read vmr profiles.
      write(*,*) version
      if (iargc() == 0) then
         write(*,*)'Enter frequency range (Fmin, Fmax)'
         read(*,*) fmin,fmax
      elseif (iargc() == 1) then
         call getarg(1, menuinputfile)
         open(10, file=menuinputfile, status='old')
         read(10,*) fmin,fmax
      else
         stop 'Usage: $gggpath/bin/mw inputfile containing selections'
      endif
      nmp=(fmax-fmin)/dnu
      if(nmp.gt.mmp) then
        write(*,*)' Increase parameter MP to ',nmp
        close(10)
        stop
      endif
      if (iargc() == 0) then
         write(*,*)'Enter target gas id & isotope (-ve = all)'
         read(*,*)kgas,kiso
         write(*,*)'Enter minimum altitude (km) along ray path '
         read(*,*)zmin
         write(*,*)'Enter zenith angle (deg) at observer '
         read(*,*)theta
      elseif (iargc() == 1) then
         read(10,*)kgas,kiso
         read(10,*)zmin
         read(10,*)theta
      endif
      close(10)
      if(theta.gt.90) then
        zobs=(zmin+dble(roc))/dsin(theta*d2r)-dble(roc)
      else
        zobs=zmin
      endif
c
      open(lunr,file='/home/toon/ggg/vmrs/gnd/gnd_summer.vmr',
     & status='old')
      read(lunr,*)nlhead,ncol
      ngas=ncol-1
      write(*,*)ngas
      if(ngas.gt.mgas) stop 'NSPECI > MSPECI'
      do j=2,nlhead
         read(lunr,'(a)') string
      end do
      do ilev=1,mlev
         read(lunr,*,end=88) z(ilev),(vmr(igas,ilev),igas=1,ngas)
         call us_std_atm(z(ilev),t(ilev),p(ilev),d(ilev))
      end do
      read(lunr,*,end=88)
      stop 'Increase parameter MLEV'
88    close(lunr)
      nlev=ilev-1
c
c   Compute Slant Columns 
      call raypath(roc,zobs,theta,z,sp,nlev)

c   Compute average (absorber-weighted) pressure and temperature
      do igas=1,ngas
         pbar(igas)=0.0
         tbar(igas)=0.0
         slcol(igas)=0.0
         sxtot(igas)=0.0
         wbar(igas)=0.0
      end do
      do ilev=1,nlev
c         write(*,*)ilev,z(ilev),sp(ilev)
         slp=d(ilev)*sp(ilev)
         do igas=1,ngas
            slpv=slp*vmr(igas,ilev)
            slcol(igas)=slcol(igas)+slpv
            pbar(igas)=pbar(igas)+p(ilev)*slpv
            tbar(igas)=tbar(igas)+t(ilev)*slpv
         end do
      end do
      do igas=1,ngas
         pbar(igas)=pbar(igas)/slcol(igas)
         tbar(igas)=tbar(igas)/slcol(igas)
         tfac(igas)=1.4388*(1/296.-1/tbar(igas))
c         write(*,*)igas,tbar(igas),pbar(igas),slcol(igas)
      end do
c
c      write(*,*)' Computing slant column abundances for each gas...'
c
c  Read through linelist, calculating absorptions for each gas
      llformat='(i2,i1,f12.6,e10.3,10x,f5.0,f5.4,f10.4,f4.2,f8.6,a24)'
      reclen=101
      llname='/home/toon/ggg/linelist/atm.161'
      do linelist=1,2
        fsib=file_size_in_bytes(lunr,llname)
        nlines=fsib/reclen
        write(*,*)' Reading '//llname, fsib
        open(lunr,file=llname,access='direct',form='formatted',
     &  status='old',recl=reclen)
        kline1=posnall(lunr,fmin,nlines)  ! index of the last line with v < NU1
        kline2=posnall(lunr,fmax,nlines)  ! index of the last line with v < NU2
c        write(*,*)kline1,kline2,nlines
        do jline=kline1+1,kline2
           read(lunr,llformat,rec=jline)
     &     molno,iso,freq,stren,abhw,sbhw,eprime,tdabhw,pshift,rotate
           sx=stren*exp((eprime-300)*tfac(molno))*slcol(molno)/dnu
           findex=(freq-fmin)/dnu  ! move line to nearest grid point
           kk=max0(1,nint(findex))
c           write(*,*)freq,index,sx
           sxtot(molno)=sxtot(molno)+sx
           wbar(molno)=wbar(molno)+abhw*sx
c           write(*,llformat)molno,iso,freq,stren
           if( molno.eq.kgas .and. (iso.eq.kiso .or. kiso.le.0) ) then
              target(kk)=target(kk)-sx
              gse(kk)=gse(kk)-sx*eprime
           else
              vac(kk,molno)=vac(kk,molno)-sx
           endif
         end do   !  jline=kline1+1,kline2
         close(lunr)
         llname='/home/toon/ggg/linelist/gct.101'
      end do   !  linelist=1,2
c
      do ip=1,nmp
         gse(ip)=gse(ip)/target(ip)
      end do
c  Convolve VAC with pressure-broadening (assumed to vary only with the Pbar)
c  Find the smallest power of 2 to accomodate NMP.
      nfft=1
      do while (nfft.lt.nmp)
         nfft=nfft+nfft
      end do
      if (nfft.gt.mmp) stop 'fringes: Increase parameter MMP'
      write(*,*)nmp,nfft
c
c  Convolve non-target lines with appropriate Lorentzians 
c  Sum them all together into array other
      do igas=1,mgas
         if(sxtot(igas).gt.0) then
            call ffak(vac(1,igas),nfft) ! Fast Fourier Analysis
            wbar(igas)=wbar(igas)/sxtot(igas)
            zz=-2*dpi*wbar(igas)*pbar(igas)/nfft/dnu
            do i=4,nfft,2
               xx=exp(i*zz)
c               vac(i-1,igas)=vac(i-1,igas)*xx  ! Real
c               vac(i,igas)=vac(i,igas)*xx      ! Imaginary
               tother(i-1)=tother(i-1)+vac(i-1,igas)*xx ! Real
               tother(i)=tother(i)+vac(i,igas)*xx       ! Imag
            end do
c            vac(2,igas)=vac(2,igas)*exp(nfft*zz)  ! Nyquist frequency
            tother(2)=tother(2)+vac(2,igas)*exp(nfft*zz) ! Nyquist
            tother(1)=tother(1)+vac(1,igas)*1            ! DC
         endif
      end do  !  igas=1,mgas
      call ffsk(tother,nfft)  ! Fast Fourier Synthesis
c
c  Convolve target lines with appropriate Lorentzian 
      call ffak(target,nfft) ! Fast Fourier Analysis
      zz=-2*dpi*wbar(kgas)*pbar(kgas)/nfft/dnu
      do i=4,nfft,2
         xx=exp(i*zz)
         target(i-1)=target(i-1)*xx  ! Real
         target(i)=target(i)*xx      ! Imaginary
      end do
      target(2)=target(2)*exp(nfft*zz)  ! Nyquist frequency
      call ffsk(target,nfft)  ! Fast Fourier Synthesis
c
c Compute solar spectrum and add to TOTHER.
         frac=0.0
         call solar_pseudo_trans_spec
     &  (19,'/home/toon/ggg/linelist/solar_dc.101',
     &  fmin-dnu,dnu,frac,spts,nmp)

c  Compute transmittances
      do ip=1,nmp
         target(ip)=exp(target(ip))
         tother(ip)=spts(ip)*exp(tother(ip))
      end do
c
c  Evaluate signal to noise ratio
      smax=0.0
      do ip=2,nmp-1
         sig_other=3-tother(ip)-tother(ip-1)-tother(ip+1)
         curv_target=abs(2*target(ip)-target(ip-1)-target(ip+1))
         curv_other=abs(2*tother(ip)-tother(ip-1)-tother(ip+1))
         snr=tother(ip)*(curv_target+1-target(ip))/
     &   (noise+curv_other+sig_other)
         if(snr .gt. smax) then
            smax=snr
            lmax=ip
         endif
      end do
      write(*,*)lmax,smax
c
c  Print out details of promising windows
      open(lunw,file='mw.out',status='unknown')
      open(lunt,file='mw.plt',status='unknown')
      write(lunw,*) 6,5
      write(lunt,*) 6,3
      write(lunw,*) version
      write(lunt,*) version
      write(lunw,'(a,3f10.2)')' Fmin, Fmax, Fdel = ',fmin,fmax,dnu
      write(lunt,'(a,3f10.2)')' Fmin, Fmax, Fdel = ',fmin,fmax,dnu
      write(lunw,*)' Gas and Isotope # = ',kgas,kiso
      write(lunt,*)' Gas and Isotope # = ',kgas,kiso
      write(lunw,*)' Minimum Altitude, Zenith Angle, Slant Column: ',
     & zmin,theta,slcol(kgas)
      write(lunt,*)' Minimum Altitude, Zenith Angle, Slant Column: ',
     & zmin,theta,slcol(kgas)
      write(lunw,*)' freq       T_target      T_other      snr       E"'
      write(lunt,*)' freq  T_target  T_other'
      do ip=2,nmp-1
         sig_other=3-tother(ip)-tother(ip-1)-tother(ip+1)
         curv_target=abs(2*target(ip)-target(ip-1)-target(ip+1))
         curv_other=abs(2*tother(ip)-tother(ip-1)-tother(ip+1))
         snr=tother(ip)*(curv_target+1-target(ip))/
     &   (noise+curv_other+sig_other)
         write(lunt,*)fmin+dnu*ip,target(ip),tother(ip)
         if(snr.gt.smax/100.and.target(ip).lt.1
     &   .and. target(ip).gt.0.5 .and. gse(ip).gt.0)
     &   write(lunw,'(f9.3,4f12.4)')
     &   ip*dnu+fmin,target(ip),tother(ip),snr,gse(ip)
      end do
      close(lunw)
      close(lunt)
      stop
      end

      subroutine raypath(roc,zobs,theta,z,sp,nlev)
c  Calculates the atmospheric ray paths associated with the levels Z(NLAY)
c
c  INPUTS
c     ROC      R*4   Radius of Curvature (of atmosphere) in km.
c     ZOBS     R*4   Observation altitude (km)
c     THETA    R*4   Solar zenith angle (deg) at observer
c     Z(NLAY)  R*4   Vector of atmospheric level altitudes (km)
c     NLAY     I*4   Number of atmospheric levels
c
c  OUTPUTS:
c     SP(NLAY) R*4   Path distances (km) associated with each level
c    
c  Notes:
c 1) Assumes that elements of Z increase monotonically
c 2) Does not include effect of refraction.

      integer*4 nlev,ilev
      real*8 d2r,robs,rmin,rup,rwas,rt,ss,swas,factor
      real*4 z(nlev),sp(nlev),roc,zobs,theta,dpi
      parameter(dpi=3.14159265359d0)
      parameter (d2r=dpi/180.0d0)

      factor=1.0d0
      robs=dble(roc)+dble(zobs)
      rt=robs*dsin(theta*d2r)  ! tangent radius
      if(theta.ge.90.0) then
         rmin=rt
         swas=0.0d0
      else
         rmin=robs
         swas=dsqrt(rmin**2-rt**2)
      endif
      write(*,*)zobs,theta,roc,robs,rt,rmin
      rwas=rmin
      do ilev=1,nlev-1
         sp(ilev)=0.0
         rup=dble(roc)+(z(ilev)+z(ilev+1))/2
         if(rup.le.rmin) then
            rwas=rup
         else
            if(rwas.lt.robs) then
               rwas=dmin1(rup,robs)
               ss=dsqrt(rwas**2-rt**2)
               factor=2.0
               sp(ilev)=factor*(ss-swas)
               swas=ss
            endif
            if(rup.gt.robs) then
               factor=1.0
               rwas=rup
               ss=dsqrt(rup**2-rt**2)
               sp(ilev)=sp(ilev)+factor*(ss-swas)
               swas=ss
            endif
         endif     !  rup.le.rmin
      end do    !  ilev=1,nlev-1
      ss=dsqrt((roc+z(nlev))**2-rt**2)
      sp(nlev)=factor*(ss-swas)
      return
      end


      subroutine us_std_atm(z,t,p,d)
c
c  Computes the 1976 US Standard Atmosphere (temperature, pressure
c  and density) at a user-supplied geometric altitude.
c
c  INPUTS:
c    Z  R*4  Geometric altitude (in km)
c
c  OUTPUTS:
c    T  R*4  Temperature (K)
c    P  R*4  Pressure (atm)
c    D  R*4  Density (molecules cm-2 km-1)
c
c  Divides the atmosphere into 8 altitude zones
c
c    Zone    H0       T0       P0         LR
c      1     0.0+    288.15 1013.25     -6.5
c      2    11.0+    216.65  226.32      0.0
c      3    20.0+    216.65   54.748     1.0
c      4    32.0+    228.65    8.8906    2.8
c      5    47.0+    270.65    1.1090    0.0
c      6    51.0+    270.65    0.66938  -2.8
c      7    71.0+    214.65    0.039564 -2.0
c      8    84.5+    186.95    0.0039814 2.0
c
      integer*4 izone,nzone
      parameter (nzone=8)
      real*4 z,h,t,p,d,mmw,gs,re,gas,avogadro,con,
     & h0(nzone),t0(nzone),p0(nzone),lr(nzone)
      parameter (gas=8.31432,gs=9.80665,avogadro=6.02217e+23,
     & mmw=28.964,re=6356.766)
      data h0/0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.5/
      data t0/288.15,216.15,216.15,228.65,270.65,270.65,214.65,186.95/
      data p0/1013.25,226.32,54.748,8.8906,1.109,.66938,.039564,.00398/
      data lr/-6.5, 0.0, 1.0, 2.8, 0.0, -2.8, -2.0, 2.0/
c
      con=mmw*gs/gas
      h=z/(1+z/re) ! convert from geometric to geopotential altitude
      do izone=nzone,2,-1
        if( h. ge. h0(izone)) go to 77
      end do
77    t=t0(izone)+lr(izone)*(h-h0(izone))
      if( lr(izone) .eq. 0.0 ) then
         p=p0(izone)*exp((h0(izone)-h)*con/t)
      else
         p=p0(izone)*(t/t0(izone))**(-con/lr(izone))
      endif
      p=100*p              ! convert from mbar to N.m-2
      d=avogadro*p/t/gas   ! molec.m-3
      p=p/101325           ! convert pressure to atm
      d=d/10               ! convert to molec.cm-2.km-1
      return
      end
