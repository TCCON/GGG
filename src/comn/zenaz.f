      subroutine zenaz(object,oblat,oblon,obalt,
     &  iy,im,id,frd,asza,azim,eorv,ervc,tplat,tplon,tpalt)
c  Computes solar zenith and azimuth angles, and the tangent point location.
c  Also computes Doppler shift components along Earth-Sun vector due to:
c  (i) the Earth's rotation, and (ii) changes in the Earth-Sun distance
c
c  Inputs:
c     object   I*4  Astronomical object (Moon=1; Sun=2)
c     oblat    R*8  Latitude of Observer (deg.)
c     oblon    R*8  Latitude of Observer (deg.)
c     obalt    R*8  Altitude of Observer (km)
c     iy       I*4  Year (e.g. 1998)
c     im       I*4  Month
c     id       I*4  Date
c     frd      R*8  Time of day (fractional days)
c
c  Outputs:
c     asza     R*8  Solar zenith angle   (deg.)
c     azim     R*8  Solar azimuth angle  (deg.)
c     eorv     R*8  Earth-Object Radial velocity (m/s)
c     ervc     R*8  Earth Rotational Velocity Component (m/s)
c     tplat    R*8  Tangent Point Latitude (deg.)
c     tplon    R*8  Tangent Point Longitude (deg.)
c     tpalt    R*8  Tangent Point Altitude (km)

      implicit none
      integer*4  object,iy,im,id

      real*8
     & sxx,cxx,d2r,
     & oblat, coblat, soblat,   !  Observer Latitude (degs.)
     & oblon,                   !  Observer Longitude (degs.)
     & obalt,                   !  Observer Altitude (km)
     & sslat, csslat, ssslat,   !  Sub-Solar Latitude (rads.)
     & sslon,                   !  Sub-Solar Longitude (rads.)
     & ddlon, sddlon, cddlon,   !  subsolar-observer longitudes
     & tplat, tplon,  tpalt,    !  Tangent Point location (degs & km)
     & asza,  casza,  sasza,    !  Astronomical Solar Zenith Angle
     & azim,   !  Solar Azimuth Angle (degs.)
     & gcrad,  !  Geocentric Radius of the Earth at OBLAT.
     & frd,    !  Time of day (fractional days)
     & pllx(2),!  Parallax: Earth-Radius/Earth-Object distance
     & ervc,   !  Earth Rotational Velocity Component (m/s)
     & eorv    !  Earth-Object Radial Velocity (m/s)
      parameter (d2r=3.14159265d0/180.0d0)
      data pllx/.0166,.000043/  ! Earth Radius / Earth-Object distance
c
c      write(*,*) 'in zenaz',object,oblat,oblon,obalt,iy,im,id,frd
c
c---------------------------------------------------------------------
c  Compute sub-solar (or sub-lunar) latitude and longitude
c  and Earth-Object Radial Velocity (eorv).
      call subsolar(object,iy,im,id,frd,sslat,sslon,eorv)
c---------------------------------------------------------------------
c  Compute astronomical solar zenith angle (ASZA).
      ssslat=dsin(sslat)
      csslat=dcos(sslat)
      soblat=dsin(oblat*d2r)
      coblat=dcos(oblat*d2r)
      ddlon=(sslon-oblon*d2r)
      sddlon=dsin(ddlon)
      cddlon=dcos(ddlon)
      casza=cddlon*csslat*coblat+ssslat*soblat
      sasza=dsqrt(1-casza**2)
c--------------------------------------------------------------------
c  Correct ASZA for parallax (due to finite Earth radius).
c     PLLX is Earth radius / Earth-Object distance
c     Parallax correction adds up to 0.002 deg in the case of the sun
c     and up to 0.9 deg in the case of the moon.
      asza=(dacos((casza-pllx(object))/
     & sqrt(1-2*pllx(object)*casza+pllx(object)**2)))/d2r
c--------------------------------------------------------------------
c  Compute azimuth angle (AZIM).
      if(sasza.gt.0.0001) then
        azim=(dacos((coblat*ssslat-soblat*csslat*cddlon)/sasza))/d2r
        if(sddlon.lt.0.0d0) azim=360.d0-azim
      else
        azim=-999.d0 !    Sun overhead - Azimuth angle ill-defined
      endif
c  Although the expression above for the solar azimuth angle can be
c  simplified to  AZIM=DACOS[(ssslat-soblat*casza)/(sasza*coblat)],
c  this formula introduces a singularity when the observer is at the
c  poles where coblat=0
c--------------------------------------------------------------------
c  Calculate the sunward Doppler shift due to the Earth's rotation
c     6378.178  is the Earth's Equatorial Radius in km.
c     3.36d-3 is the Earth's oblateness
c     7.292116d-5 is the Earth's rotational velocity (radians/s)
      gcrad=1000.d0*(obalt+6378.178d0*(1.d0-3.36d-3*soblat**2)) ! Geocentric Radius
      ervc=-(7.2722d-5*gcrad*coblat*csslat*sddlon)       ! velocity toward sun
c     -ve values imply distance to sun is getting smaller
c-----------------------------------------------------------------------
c  Calculate the latitude and longitude of the unrefracted tangent point.
c  If the solar zenith angle at the observer < 90 deg. then a virtual
c  tangent point exists behind the observer. Note that the location (lat,long)
c  of the tangent point is independent of the altitude of the observer.
c  Also note that the (unrefracted) distance to the tangent point is simply
c                -(GCRAD+OBALT)*COS(ASZA)
c  which is -ve in the case of a virtual tangent point.
      sxx=coblat*dsin(oblon*d2r)-casza*csslat*dsin(sslon)
      cxx=coblat*dcos(oblon*d2r)-casza*csslat*dcos(sslon)
      tplat=dasin((soblat-casza*ssslat)/sasza)/d2r
      tplon=datan2(sxx,cxx)/d2r
      tpalt=sasza*(obalt+gcrad/1000)-gcrad/1000
c-----------------------------------------------------------------------
      return
      end
