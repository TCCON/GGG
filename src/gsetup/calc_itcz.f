      subroutine calc_itcz(alon_obs,idoy_obs,itcz_lat,itcz_width)
c  Calculate the ITCZ latitude based on the day of year and longitude.

c  Inputs:
c         alon_obs        R*8   Longitude of site/observation (deg)
c         idoy_obs        I*4   Day of year of observation
c
c Outputs:
c         itcz_lat        R*8   Latitude of the centre of the itcz (deg) 
c         itcz_width      R*8   Latitude width of the itcz  (deg)
c
c the outputs will be applied in the resample_vmrs_at_effective_altitudes
c 
c written by minqiang.zhou@aeronomie.be   2017-01-30
c
c Modified by GCT to perform longitude interpolation, use DOY
c (instead of month), and tabulte onto a 15 deg longitude grid.

      implicit none
      integer*4 nlon,ilon,idoy_obs
      parameter (nlon=24)  ! 24 x 15 deg = 360 deg
      real*8 vlat_jul(0:nlon),vlat_jan(0:nlon),vwidth(0:nlon),
     &alon_obs,xlon,fr,janlat,jullat,itcz_lat,itcz_width,pi
c==================================================================
c  Latitudes and widths of ITCZ in July and January (every 15 deg from 0 to 360)   
c  First and last points (O and 360 deg) are repeated to simplify interpolation.
c  from wiki https://en.wikipedia.org/wiki/Intertropical_Convergence_Zone#/media/File:ITCZ_january-july.png  
c  Assume peak Northward excursion of ITCZ is mid-July (idoy=198)
c  Assume peak Southward excursion of ITCZ is mid-Jan (idoy=15)
      vlat_jul=(/+16,+20,+24,+27,+29,+30,+30,+30,+29,+26,+21,+20,
     &           +19,+19,+19,+19,+18,+16,+11, +5, +1, +4, +8,+12,+16/) 
      vlat_jan=(/ +3, +5, -8,-14,-12, -8, -5, -4, -5, -7,-10,-13,
     &           -14,-14,-10, -4, -0, +1, +0, -2, -6, -8, -6, -2, +3/)
      vwidth=  (/+11,+10, +9, +8, +8, +9,+10,+11,+12,+12,+11,+10,
     &            +9, +8, +7, +6, +7, +8,+11,+13,+14,+14,+13,+12,+11/)

      pi=2*dacos(0.0d0)

c Longitude interpolation of the ITCZ latitudes & widths.
      if(alon_obs.lt.0) alon_obs=alon_obs+360
      xlon=nlon*alon_obs/360.0
      ilon=int(xlon)
      fr=xlon-ilon
      ilon=mod(ilon,nlon)
      
      itcz_width=(1-fr)*vwidth(ilon)+fr*vwidth(ilon+1)
      janlat=(1-fr)*vlat_jan(ilon)+fr*vlat_jan(ilon+1)
      jullat=(1-fr)*vlat_jul(ilon)+fr*vlat_jul(ilon+1)
      itcz_lat = 0.5*(jullat+janlat -
     & (jullat-janlat)*dcos(2*pi*(idoy_obs-15)/365.25))

      return
      end
