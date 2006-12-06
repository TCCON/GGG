      subroutine setvmr(vmr,mgas,z,t,p,h2ovmr,nlev,iyr,iday,
     & uttime,oblat,oblon,obalt,tout,pout,hout,ptrop_atm,ppbl_atm)
c
c  Purpose:    Allows a vmr profile to be modified
c
c  Inputs:     Everything
c
c  Outputs:    vmr
c
      implicit none

      integer*4
     & mgas,       ! Maximum allowed number of gases
     & nlev, ilev, ! Number of atmospheric levels (70?)
     & iyr,        ! year (e.g. 2005)
     & iday        ! Day of year (1-366)

      real*4 vco2ref,dz2,vmr0,age,apbl,atrop,atoa,seasonal,sdma,
     & vmr(mgas,nlev),   ! buffer for vmr's
     & p(nlev),          ! pressures of levels (atm.)
     & t(nlev),          ! temperatures of levels (K)
     & z(nlev),          ! altitudes of levels (km)
     & h2ovmr(nlev)      ! H2O vmr profile from .mod file

      real*8 pi,
     & uttime,           ! UT Time (hours)
     & oblat,            ! Observation Latitude (deg)
     & oblon,            ! Observation Longitude (deg)
     & obalt,            ! Observation Altitude (km)
     & tout,             ! Outside temperature (c)
     & pout,             ! Outside pressure (mbar)
     & hout,             ! Outside RH (%)
     & ptrop_atm,        ! Tropopause pressure (atm)
     & ppbl_atm          ! PBL pressure (atm)

      pi=2*dacos(0.0d0)
      vco2ref=0.000380

      if(ptrop_atm.gt.0) then  ! fudge CO2 only if ptrop is non-zero.

c  If no PBL pressure-altitude supplied, assume that the daytime high-latitude
c  PBL varies from 850 mbar in winter to 650 mbar in summer.
         if( ppbl_atm.le.0.0) ppbl_atm=(0.75-0.10*sin(2*pi*oblat/360)*
     &    sin(2*pi*(iday+237)/365.25))
c      write(*,*)'setenv:',ptrop_atm,ppbl_atm
c
c  Calculate age of air for each atmospheric level.
         apbl=0.2  ! Age of air at the PBL (years)
         atrop=0.4 ! Age of air at the tropopause (years)
         atoa=5.5  ! Age of air in the upper strat (years)
         do ilev=1,nlev
            if(p(ilev).gt.ppbl_atm) then
               age=apbl*(p(1)-p(ilev))/(p(1)-ppbl_atm)
            elseif(p(ilev).gt.ptrop_atm) then
               age=apbl+(atrop-apbl)*(ppbl_atm-p(ilev))/
     &         (ppbl_atm-ptrop_atm)
            else
               age=atoa-(atoa-atrop)*sqrt(p(ilev)/ptrop_atm)
            endif
c            write(*,*)ilev,p(ilev),age
            vmr0=vco2ref*(1-0.005*age) ! assume CO2 increases at 0.5%/year
            seasonal=0.01*sin(2.0*oblat*(1-oblat/720)*pi/180)*  ! Amplitude of seasonal variation
     &      exp(oblat/45)
            sdma=sin(2*pi*(float(iday+75)/365.25-age))     ! Sinusoidal seasonal variation
            sdma=(1.45-exp(-(1.05*sdma)))                  ! Non-sinusoidal seasonal cycle
            vmr(2,ilev)=vmr0*(1+seasonal*exp(-(age/0.25))*sdma)
         end do

      endif  ! if(ptrop_atm.gt.0.0)
c
c  Factor 0.01*sin(d2r*2.0*oblat*(1-oblat/720))*exp(oblat/45) estimates amplitude of
c  seasonal cycle. It peaks at +3.7% at 70N, is zero at equator, and
c  peaks again at -0.45% at 30S.
c
c  Factor exp(-age/0.25) accounts for decay in seasonal cycle with age of air.
c
c  Factor (1.45-exp(-1.05*sdma)) accounts for fact that seasonal variation
c  of CO2 is not sinusoidal. The drawdown happens very quickly over about
c  8 weeks in early summer. This modifies the sinudoidal variation of sdma.


c     Substitution of model H2O profile
      if(h2ovmr(1).gt.0.0) then ! An h2o profile was found in the .mod file
         if ( pout.gt.(0.6) ) then ! ground-based observations
            do ilev=1,nlev
               vmr(1,ilev)=h2ovmr(ilev)  ! H2O
               vmr(49,ilev)=h2ovmr(ilev)*0.16*(8.0+log(h2ovmr(ilev)))  ! HDO
            end do
         endif
      endif
c
c  The factor 0.16*(8.0+log(h2ovmr(ilev))) represents the isotopic fractionation
c  of HDO/H2O.

      return
      end
