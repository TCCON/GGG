      real*4 function compute_ztrop_gct(nlev,z,t)
c
c  Finds/computes the first/lowest altitude at which dT/dz > -2 K/km
c  having previously found a lower altitude at which dT/dz < -2 K/km
c  
c  Smoothes raw data to ~2 km vertical resolution, avoiding the need
c  to check whether dT/dz falls back below -2 k/km within 2km of
c  current level.
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
c  Bins the raw Z/T data onto a uniform verical grid of spacing dz.
c  This has the benefit of averaging and smoothing the raw data,
c  which reduces the liklihood of a false tropopsuse being detected
c  in oversampled or noisy data.  This also means that it doesn't
c  matter whether the raw data is sorted by altitude or not.
c
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
c  having previously been less than -2
c   -2 = LRm + [LRp-LRm][zt-zm]/[zp-zm]
c
c  Hence
c   zt = zm + [zp-zm][-2-LRm]/[LRp-LRm]

c  Implementation Notes:
c  1) lrz = Lapse Rate is defined here as dT/dz, whereas elsewhere in
c   the meteorological literature it is typically defined at -dT/dz
c  2) qrz = Lapse rate per level (qrz=lrz when dz=1.0)
c  2) Assumes that lapse rate is a linear & continuous function of Z
c  3) Can never find a tropopause below zlo or above zhi
c  4) Input vectors unaltered
c  5) compute_ztrop = -99.9 if no tropopause found
c    confusing an inversion at top of PBL for the tropopause.

c  Lapse rates are computed at the altitudes mid-way between levels
c  by LRZ = [t(i+1)-t(i)]/[z(i+1)-z(i)]
c  at z = zlr=[z(i+1)+z(i)]/2

      implicit none
      integer*4 nlev,ilev,klev,iflag,
     & kmin,kmax,kk,mk
      parameter (mk=120)
      real*4 z(nlev),t(nlev),tt(mk),tw(mk),wt,
     & zlo,zhi,
     & dz,uu,
     & qrz,qrzwas,lrthresh,qrthresh
      parameter(zlo=3.5,zhi=25.0,
     & lrthresh=-2.0,  ! Lapse Rate Threshold (K/km)
     & dz=1.2,
     & qrthresh=lrthresh*dz)
   
      do klev=1,mk
         tw(klev)=0.0
         tt(klev)=0.0
      end do

c  Apportion each temperature measurement between two
c  adjacent levels. Also find min/max levels needed.
      kmin=mk
      kmax=1
      do ilev=1,nlev
         if( z(ilev).ge.zlo .and. z(ilev).lt.zhi ) then
            uu=z(ilev)/dz
            kk=1+int(uu)
            wt=kk-uu
            tw(kk)=tw(kk)+wt
            tw(kk+1)=tw(kk+1)+1-wt
            tt(kk)=tt(kk)+wt*t(ilev)
            tt(kk+1)=tt(kk+1)+(1-wt)*t(ilev)
            if(z(ilev).lt.dz*(kmin-1)) kmin=1+int(z(ilev)/dz)
            if(z(ilev).gt.dz*(kmax-1)) kmax=1+int(z(ilev)/dz)
         endif
      end do

      iflag = 0   !  Flag indicating whether LR < LRT has been found yet
c   IFLAG=1 means that dT/dz < -2 K/km has been encountered above 600 mbar
      do klev=kmin,kmax         ! loop over levels
c         write(*,*) 2,klev,dz*(klev-1),iflag,tt(klev)/tw(klev)
         if(tw(klev+1).le.0.0) stop  'Hole encountered: increase dz'
         qrz=tt(klev+1)/tw(klev+1)-tt(klev)/tw(klev)    ! Lapse Rate at ZLR
         if( iflag.gt.0 .and. qrz.gt.qrthresh ) go to 99  ! Success
         qrzwas=qrz
         if(qrz.le.qrthresh) iflag=1
      end do
c      write(*,*) 'Warning: Tropopause not found'
      compute_ztrop_gct=-99.9
      return
c   zt = zm + [zp-zm][LRt-LRm]/[LRp-LRm]
c   zt = zm + dz.[LRt-LRm]/[LRp-LRm]
99    compute_ztrop_gct=dz*(klev-1.5+(qrthresh-qrzwas)/(qrz-qrzwas))
c      write(*,*)'compute_ztrop_gct=',compute_ztrop_gct
      return
      end
