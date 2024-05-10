      real*4 function gravity(gdlat_deg,altit_km)
c
c  Input Parameters:
c      gdlat_deg   r*4  Geodetic Latitude (deg)
c      altit_km    r*4  Geometric Altitude (km)
c
c  Output Parameter:
c      gravity     r*4  Effective Gravitational Acceleration (m/s2)
c
c  Computes the effective Earth gravity at a given latitude and altitude.
c  This is the sum of the gravitational and centripital accelerations.
c  These are based on equation I.2.4-(17) in US Standard Atmosphere 1962
c  The Earth is assumed to be an oblate ellipsoid, with a ratio of the
c  major to minor axes = sqrt(1+cc) where cc=.006738
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

      real*4  d2r,gm,omega,cc,shc,eqrad,gdlat_deg,gdl_rad,altit_km,
     & gclat,         ! geocentric latitude (radians)
     & sgcl,cgcl,
     & radius,        ! radial distance (metres)
     & ff,hh,ge       ! scratch variables

      parameter(
     & d2r=3.14159265/180,! Conversion from degrees to radians
     & gm=3.9862216e+14,  ! Gravitational constant times Earth's Mass (m3/s2)
     & omega=7.292116E-05,! Earth's angular rotational velocity (radians/s)
     & cc=0.006738,       ! (a/b)**2-1 where a & b are equatorial & polar radii
     & shc=1.6235e-03,    ! 2nd harmonic coefficient of Earth's gravity field 
     & eqrad=6378178.)    ! Equatorial Radius (m)
c
c  Convert from geodetic latitude (GDLAT) to geocentric latitude (GCLAT).
      gdl_rad=d2r*gdlat_deg
      gclat=atan(tan(gdl_rad)/(1+cc))  ! radians
c  On computers which crash at the poles using the equation above, try
c      gclat=gdl_rad-cc*sin(gdl_rad)*cos(gdl_rad)/(1+cc*cos(gdl_rad)**2)

      sgcl=sin(gclat)
      cgcl=cos(gclat)
      radius=1000*altit_km + eqrad/sqrt(1.+cc*sgcl**2) ! meters
      ff=(radius/eqrad)**2
      hh=radius*omega**2
      ge=gm/eqrad**2                       ! = gravity at Re
      gravity=(ge*(1-shc*(3*sgcl**2-1)/ff)/ff-hh*cgcl**2)*
     & (1+0.5*(sgcl*cgcl*(hh/ge+2*shc/ff**2))**2)
      return
      end
