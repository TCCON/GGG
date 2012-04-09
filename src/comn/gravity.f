      real*4 function gravity(gdlat,altit)
c
c  Input Parameters:
c      gdlat       r*4  GeoDetric Latitude (degrees)
c      altit       r*4  Geometric Altitude (km)
c
c  Output Parameter:
c      gravity     r*4  Effective Gravitational Acceleration (m/s2)
c
c  Computes the effective Earth gravity at a given latitude and altitude.
c  This is the sum of the gravitational and centripital accelerations.
c  These are based on equation I.2.4-(17) in US Standard Atmosphere 1962
c  The Earth is assumed to be an oblate ellipsoid, with a ratio of the
c  major to minor axes = sqrt(1+con) where con=.006738
c  This eccentricity makes the Earth's gravititational field smaller at
c  the poles and larger at the equator than if the Earth were a sphere
c  of the same mass. [At the equator, more of the mass is directly
c  below, whereas at the poles more is off to the sides). This effect
c  also makes the local mid-latitude gravity field not point towards
c  the center of mass.
c
c  The equation used in this subroutine agrees with the International
c  Gravitational Formula of 1967 (Helmert's equation) within 0.005%.
c
c  Interestingly, since the centripital effect of the Earth's rotation
c  (-ve at equator, 0 at poles) has almost the opposite shape to the
c  second order gravitational field (+ve at equator, -ve at poles),
c  their sum is almost constant so that the surface gravity could be
c  approximated (.07%) by the simple expression g=0.99746*GM/radius^2,
c  the latitude variation coming entirely from the variation of surface
c  r with latitude. This simple equation is not used in this subroutine.

      implicit none
      include "../ggg_const_params.f"

      real*4  d2r,gm,omega,con,shc,eqrad,gdlat,altit,
     & gclat,         ! geocentric latitude (radians)
     & radius,        ! radial distance (metres)
     & ff,hh,ge       ! scratch variables
      parameter(
     & d2r=spi/180,       ! Conversion from degrees to radians
     & gm=3.9862216e+14,  ! Gravitational constant times Earth's Mass (m3/s2)
     & omega=7.292116E-05,! Earth's angular rotational velocity (radians/s)
     & con=.006738,       ! (a/b)**2-1 where a & b are equatorial & polar radii
     & shc=1.6235e-03,    ! 2nd harmonic coefficient of Earth's gravity field 
     & eqrad=6378178.)    ! Equatorial Radius (m)
c
c  Convert from geodetic latitude (GDLAT) to geocentric latitude (GCLAT).
      gclat=atan(tan(d2r*gdlat)/(1+con))  ! radians
c  On computers which crash at the poles try the following expression
c  gclat=d2r*gdlat-con*sin(d2r*gdlat)*cos(d2r*gdlat)/(1+con*cos(d2r*gdlat)**2)
      radius=1000*altit+eqrad/sqrt(1.+con*sin(gclat)**2)
      ff=(radius/eqrad)**2
      hh=radius*omega**2
      ge=gm/eqrad**2                       ! = gravity at Re
      gravity=(ge*(1-shc*(3*sin(gclat)**2-1)/ff)/ff-hh*cos(gclat)**2)
     & *(1+0.5*(sin(gclat)*cos(gclat)*(hh/ge+2*shc/ff**2))**2)
      return
      end
