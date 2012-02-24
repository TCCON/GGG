      function compute_seasonal_cycle(kgas,zobs,ztrop,alat_obs,fryr)
c
c  Purpose: 
c         Modifies the a priori vmr profiles on a gas-by-gas basis
c         to account for the season of the observation/model.
c
c  Inputs:
c         kgas,             I*4   Gas ID (1=H2O; 2=CO2; 3=O3; etc)
c         zobs,             R*8   Observation/Site altitude
c         ztrop,            R*8   Observation/Site tropopause altitude
c         alat_obs          R*8   Latitude of site/observation
c         fryr,             R*8   Fraction of the year
c
c  Outputs:
c         compute_seasonal_cycle    R*4   Seasonal Fractional Adjustment
c
c  Note the input VMR matrix is over-written by the output vmr matrix.
c
      implicit none
      integer*4 kgas,ngas
      parameter (ngas=6)

      real*4 calc_aoa,compute_seasonal_cycle

      real*8 seasonal(ngas),
     & fryr,ztrop,alat_obs,zobs,sv,svl,svnl,
     & sca,aoa,twopi


      twopi=4.0d0*dacos(0.0d0)

c  Amplitude at surface at 50N
      seasonal(1)=0.00   ! H2O
      seasonal(2)=0.008  ! CO2
      seasonal(3)=0.00   ! O3
      seasonal(4)=0.00   ! N2O
      seasonal(5)=0.25   ! CO
      seasonal(6)=0.00   ! CH4

      compute_seasonal_cycle=1.0
      if(kgas.gt.ngas) return
      aoa=calc_aoa(alat_obs,zobs,ztrop)
c      write(*,*)alat_obs,zobs,ztrop,aoa,exp(-aoa/0.30)
      if(kgas.eq.2) then
      sv=dsin(twopi*(fryr-0.84-aoa))  ! seasonal variation
      svnl=sv+1.25*exp(-((alat_obs-75)/40)**2)*(0.5-sv**2)  ! seasonal variation
      sca=svnl*exp(-aoa/0.205)*
     & (1+1.25*exp(-((alat_obs-75)/48)**2)*(zobs+6.0)/(zobs+1.2))
      else
        sv=dsin(twopi*(fryr-0.89))      ! basic seasonal variation
        svl=sv*(alat_obs/15)/sqrt(1+(alat_obs/15)**2) ! latitude dependence 
        sca=svl*exp(-aoa/1.60)              ! altitude dependence
      endif
      compute_seasonal_cycle=1.0+sca*seasonal(kgas)

c      sdma=sin(2*pi*(float(iday+75)/365.25-age))
c      sdma=(1.45-exp(-(1.11*sdma)))
c      co2vmr(ilev)=vmr0*(1+fasc*exp(-(age/0.25))*sdma)
      return
      end
