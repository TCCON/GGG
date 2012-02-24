      real*4 function calc_aoa(alat,zobs,ztrop)
c
c  Computes the Age of Air at particular location  (Altitude, Latitude)
c  relative to the surface at 50N, where Age=0
c  
c  Inputs:
c       ALAT   R*8  Latitude (deg)
c       ZOBS   R*8  Altitude (km)
c      ZTROP   R*8  Tropopause Altitude (km)
c  
c Output:
c    CALC_AOA  R*8  Calculated Age for Air


      implicit none
      real*8 alat,zobs,ztrop,aa

      aa=alat/21
      calc_aoa=0.320-0.085*exp(-((alat-46)/19)**2)
     & -0.27*exp(-1.60*zobs/(zobs+ztrop))*aa/sqrt(1+aa**2)
      if(zobs.gt.ztrop) calc_aoa=calc_aoa+5.5*(zobs-ztrop)/zobs
      return
      end
