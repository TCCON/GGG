      subroutine simulate_co2_vmr(iyr,iday,uttime,oblat,oblon,
     & ztrop,zpbl,z,co2vmr,nlev)
c
c  Purpose:   Generates an a priori CO2 vmr profile for any time/latitude
c
c  Inputs: 
c          IYR         I*4   Year
c          IDAY        I*4   Day of year
c          uttime      R*8   UT Time (hours)
c          oblat       R*8   Observation Latitude (deg)
c          ptrop_atm   R*8   Tropopause pressure (atm)
c          ppbl_atm    R*8   PBL pressure (atm)
c          z(nlev)     R*4   Altitude (km) of levels
c          p(nlev)     R*4   Pressure (atm) of levels
c          nlev        I*4   Number of levels
c
c  Outputs:   
c          co2vmr(nlev) R*4  CO2 vmr profile
c
c
c  Comments:
c    - Assumes a linear increase in CO2
c    - Does not currently do any dirurnal variations
c      Nor any longitudinal variation.
c    - Uses user-supplied values for the pressures at
c      the top of the PBL and the tropopause. If these
c      are zero, it makes up something.
c
c  0.01*sin(d2r*2.0*oblat*(1-oblat/720))*exp(oblat/45) is amplitude of
c  seasonal cycle. It peaks at +3.7% at 70N, is zero at equator, and
c  peaks again at -0.45% at 30S.
c
c  Factor exp(-age/0.25) accounts for decay in seasonal cycle with age of air.
c
c  Factor (1.45-exp(-1.05*sdma)) accounts for fact that seasonal variation
c  of CO2 is not sinusoidal. The drawdown happens very quickly over about
c  8 weeks in early summer. This modifies the sinudoidal variation of sdma.

      implicit none
      include "../ggg_const_params.f"

      integer*4
     & nlev, ilev, ! Number of atmospheric levels (70?)
     & iyr,        ! year (e.g. 2005)
     & iday        ! Day of year (1-366)

      real*4 vco2ref,vmr0,age,apbl,atrop,atoa,fasc,sdma,roi,
     & co2vmr(nlev),   ! buffer for vmr's
c     & p(nlev),        ! pressures of levels (atm.)
     & z(nlev)         ! Altitudes of levels

      real*8 
     & sslat,
     & uttime,           ! UT Time (hours)
     & oblat,            ! Observation Latitude (deg)
     & oblon,            ! Observation Longitude (deg)
     & ztrop,            ! Tropopause Altitude (km)
     & zpbl              ! PBL Altitude (km)

      roi=0.005   ! rate of increase of CO2 (=0.5%/year)
      vco2ref=0.000380*(1+roi*(iyr+iday/365.25-2005.0))

c  If no PBL altitude supplied, assume that the daytime high-latitude
c  PBL varies from 850 mbar in winter to 650 mbar in summer.
c  PBL varies from 1.5 km in winter to 3.5 km in summer.
c      if(zpbl.le.0.0) ppbl_atm=(0.70-0.15*cos(4*dpi*oblat/360)
c     & -0.10*sin(2*dpi*oblat/360)*sin(2*dpi*(iday-110)/365.25))
      if(zpbl.le.0.0) zpbl=(2.5+1.0*cos(4*dpi*oblat/360)
     & +1.0*sin(2*dpi*oblat/360)*sin(2*dpi*(iday-110)/365.25))
c
c  If no topopause pressure supplied, assume that it varies
c  from 90 mbar in the tropics to 280 mbar at high latitudes
c  from 16km in the tropics to 9km at high latitudes
      if(ztrop.le.0.0) then
          sslat=12*sin(2*dpi*(iday-120)/365.25)
c          ptrop_atm=0.28-0.19*exp(-((oblat-sslat)/35)**2)
          ztrop=9.0+7.0*exp(-((oblat-sslat)/35)**2)
      endif
c
c  Calculate age of air for each atmospheric level.
      apbl=0.2  ! Age of air at the PBL (years)
      atrop=0.4 ! Age of air at the tropopause (years)
      atoa=5.5  ! Age of air in the upper strat (years)
c
c  fasc = Fractional Amplitude of Seasonal Cycle (at surface)
      fasc=0.01*sin(2.0*oblat*(1-oblat/720)*dpi/180)*exp(oblat/45)
c
      do ilev=1,nlev
         if(z(ilev).lt.zpbl) then  ! below the PBL
            age=apbl*(z(1)-z(ilev))/(z(1)-zpbl)
         elseif(z(ilev).lt.ztrop) then  ! below the trop
            age=apbl+(atrop-apbl)*(zpbl-z(ilev))/
     &      (zpbl-ztrop)
         else                               ! in the stratosphere
c            age=atoa-(atoa-atrop)*sqrt(p(ilev)/ptrop_atm)
            age=atoa-(atoa-atrop)*exp(-(z(ilev)-ztrop)/15.)
         endif
         vmr0=vco2ref*(1-roi*age) ! assume CO2 increases at 0.5%/year
c  sdma = 
         sdma=sin(2*dpi*(float(iday+75)/365.25-age))    
         sdma=(1.45-exp(-(1.11*sdma)))       
         co2vmr(ilev)=vmr0*(1+fasc*exp(-(age/0.25))*sdma)
         if(z(ilev).gt.90) then
            co2vmr(ilev)=co2vmr(ilev)/sqrt(1+((z(ilev)-90.)/7.5)**3)
         endif
      end do
      return
      end
