      real*4 function compute_ztrop(nlev,z,t)
c
c  Finds/computes the first/lowest altitude at which dT/dz > -2 K/km
c  having previously found a lower altitude at which dT/dz < -2 K/km
c  
c  Does not bother checking if dT/dz falls back below -2 k/km within 2km
c  of current level,  since levels are sometimes > 2 km apart.
c
c   Inputs:
c        nlev    I*4   Number of levels
c       z(nlev)  R*4   vector of Altitudes
c       t(nlev)  R*4   vector of Temperatures
c
c   Output:
c       compute_ztrop  R*4  Tropopause altitude
c
c
c  Theoretical Basis
c  Computes the lapse rate (dT/dz) at the mid-way points between the
c  levels. When this exceeds the threshold (-2 K/km) it estimates the
c  altitude at which dT/dz was exactly -2 k/km, based on the assumption
c  that dT/dZ varies linearly between the mid-way points that bracket
c  the threshold value. Let:
c  LRp = (T(i+1)-T(i)]/[z(i+1)-z(i)] be the lapse rate at zp=[z(i+1)+z(i)]/2
c  LRm = (T(i)-T(i-1)]/[z(i)-z(i-1)] be the lapse rate at zm=[z(i)+z(i-1)]/2
c
c  Assume lapse rate is linear in altitude between zm and zp
c   LR(z) = LRm + [LRp-LRm][z-zm]/[zp-zm]
c
c  So the tropopause is the altitude at which LR(z) = -2
c   -2 = LRm + [LRp-LRm][zt-zm]/[zp-zm]
c
c  Hence
c   zt = zm + [zp-zm][-2-LRm]/[LRp-LRm]

c  Implementation Notes:
c  1) lrz = Lapse Rate is defined here as dT/dz, whereas elsewhere in
c   the meteorological literature it is typically defined at -dT/dz
c  2) Assumes that lapse rate is a linear & continuous function of Z
c  3) Can never find a tropopause below ~z(lmin)
c  4) Input vectors unaltered
c  5) compute_ztrop = -99.9 if no tropopause found
c  6) Doesn't start looking for tropopause until level lmin=6 to avoid
c    confusing an inversion at top of PBL for the tropopause.

c  Lapse rates are computed at the altitudes mid-way between levels
c  by LRZ = [t(i+1)-t(i)]/[z(i+1)-z(i)]
c  at z = zlr=[z(i+1)+z(i)]/2

      implicit none
      integer*4 nlev,ilev,iflag,lmin
      real*4 z(nlev),t(nlev),lrz,lrzwas,zlr,zlrzwas,lrthresh
      parameter(
     & lrthresh=-2.0,  ! Lapse Rate Threshold (K/km)
     & lmin=6)      ! Lowest level at which Ztrop can be found
   
      iflag = 0   !  Flag indicating whether LR < LRT has been found yet
c   IFLAG=0 is initial condition
c   IFLAG=1 means that dT/dz < -2 K/km has been encountered above 600 mbar

      zlrzwas=0.5*(z(lmin-2)+z(lmin-1))        
      lrzwas=(t(lmin-1)-t(lmin-2))/(z(lmin-1)-z(lmin-2))
      do ilev=lmin,nlev                ! loop over levels
         zlr=0.5*(z(ilev-1)+z(ilev))        ! Altitude at which lapse rate = lrz
         lrz=(t(ilev)-t(ilev-1))/(z(ilev)-z(ilev-1)) ! Lapse Rate at ZLR
c         write(*,*)ilev,z(ilev),t(ilev),zlr,lrz,iflag
         if( iflag.gt.0 .and. lrz.gt.lrthresh ) go to 99  ! Success
         lrzwas=lrz
         zlrzwas=zlr
         if(lrz.le.lrthresh) iflag=1
      end do
      write(*,*) 'Warning: Tropopause not found'
      compute_ztrop=-99.9
      return
c   zt = zm + [zp-zm][LRt-LRm]/[LRp-LRm]
99    compute_ztrop=zlrzwas+(zlr-zlrzwas)*(lrthresh-lrzwas)/(lrz-lrzwas)
c      write(*,*)'compute_ztrop=',compute_ztrop
      return
      end
