      subroutine apply_seasonal_cycle(nlev,z,mgas,ngas,seacycle,
     & ztrop,alat_obs,fryr,apvmr)
c
c  Purpose: 
c      Modifies the a priori vmr profiles on a gas-by-gas basis
c      to account for the season of the observation/model.
c
c  Inputs:
c      nlev,             I*4   Number of atmospheric levels
c      z(nlev),          R*4   Level altitudes
c      mgas,             I*4   Max number of gases
c      ngas,             I*4   Number of gases
c      seacycle(ngas)    R*4   Amplitude of seasonal cycles (from .vmr header)
c      ztrop,            R*8   Observation/Site tropopause altitude
c      alat_obs          R*8   Latitude of site/observation
c      fryr,             R*8   Fraction of the year
c      apvmr(mgas,nlev)  R*4   Incoming vmr profiles
c
c  Outputs:
c      apvmr(mgas,nlev)  R*4   Outgoing vmr profiles
c
c  Note the input VMR matrix is over-written by the output vmr matrix.

      implicit none
      integer*4 mgas,jgas,ngas,nlev,ilev

      real*4 seacycle(ngas),apvmr(mgas,nlev),z(nlev),
     & sv,svl,svnl,sca,twopi

      real*8 fryr,ztrop,alat_obs,aoa,calc_aoa

      twopi=4*acos(0.0)

c  Amplitude at surface at 50N
c      seacycle(1)=0.000  ! H2O
c      seacycle(2)=0.0081 ! CO2
c      seacycle(3)=0.000  ! O3
c      seacycle(4)=0.000  ! N2O
c      seacycle(5)=0.250  ! CO
c      seacycle(6)=0.004  ! CH4
c      do jgas=7,ngas
c         seacycle(jgas)=0.0
c      enddo

      do ilev=1,nlev
c         z8=dble(z(ilev))
         aoa=calc_aoa(alat_obs,dble(z(ilev)),ztrop)
         do jgas=1,ngas
            if(jgas.eq.2) then
               sv=sin(twopi*sngl(fryr-0.834-aoa))  ! seasonal variation GCT
               svnl=sv+1.80*exp(-(sngl(alat_obs-74)/41)**2)*(0.5-sv**2)
c  Note: 0.5-sin^2(x) = 0.5*cos(2x), so this is the twice-per-year term (semi-annual)
c  Faster to square the sin(x) function that we have already, than compute cos(2x)
               sca=svnl*exp(-sngl(aoa/0.20))*
     &         (1+1.33*exp(-(sngl(alat_obs-76)/48)**2)*
     &         (z(ilev)+6.)/(z(ilev)+1.4))
            else
               sv=sin(twopi*sngl(fryr-0.78))       ! basic seasonal variation KMS
c               sv=sin(twopi*(fryr-0.75))       ! basic seasonal variation GCT
               svl=sv*sngl(alat_obs/15)/sqrt(1+sngl(alat_obs/15)**2) ! latitude dependence 
               sca=svl*exp(-sngl(aoa/0.85))           ! altitude dependence KMS
c               sca=svl*exp(-sngl(aoa/0.5))           ! altitude dependence GCT
            endif
            apvmr(jgas,ilev)=apvmr(jgas,ilev)*(1+sca*seacycle(jgas))
         end do
      end do

      return
      end
