      real*8 function dsbydt(asza)
c  Estimates the fractional change of slant column for a +1 degree change
c  in zenith angle for a gas having a scale height of SCHT.
c  Can be used for refracted (-ve) or astronomical (+ve) zenith angles.
c  This can then be used to determine the zenith angle correction from a
c  .TAV file simply by setting  asza = asza + (totcon-1)/dsbydt(asza)
c
c  This approximation to dsbydt is independent of the atmospheric model
c  or the observation altitude and so is generally only accurate to ~20%.
c  When tangent altitudes fall below 14 km and astronomical zenith angles
c  are employed, refractive effects cause the actual dsbydt to be much
c  smaller than that calculated, leading to slow convergence for ASZA.
c  However, if refracted zenith angles are used, there is no problem.
c
c  Note that for ASZA << 90; dsbydt --> tan(d2r*asza)
c            For ASZA =  90; dsbydt --> 1/del
c            For ASZA >> 90; dsbydt --> R*cos(d2r*asza)/scht
c  And thus has all the correct behavior in these limits

      real*8 radius,scht,piby2,aa,asza,del,dpi
      parameter (dpi=3.14159265359d0)
      parameter (piby2=dpi/2,radius=6378,scht=7.0)
c
      del=sqrt(piby2*scht/radius)
      aa=abs(asza)*piby2/90
      if(aa.lt.piby2) then
        dsbydt=dsin(aa)/(del*exp(7.8*(aa-piby2))+dcos(aa))
      else
        dsbydt=(1-piby2*dcos(aa)/del)/del
      endif
        dsbydt=dsbydt*piby2/90
      return
      end
