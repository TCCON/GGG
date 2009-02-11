c                           xloc
c==================================================================
c
c       This subroutine returns the latitude, longitude, and altitude
c       corresponding to certain "standard" location codes read in
c       from the ftp.in file.  The standard location codes are:
c
c       MCM    McMurdo Station, Antarctica
c       JPL    Jet Propulsion Laboratory, California
c       MES    Mesa, Jet Propulsion Laboratory, California
c       TMO    Table Mountain Observatory, California
c       ARC    Ames Research Center, California
c       PAL    Palestine, Texas
c       FTS    Fort Sumner, New Mexico
c       DAG    Daggett, California
c       MTB    Mt. Barcroft Research Station, California
c       LYL    Lynn Lake, Manitoba
c       FAI    Fairbanks, Alaska
c       ESN    Esrange, Sweden
c
c       If a code other than these is passed, the value of 0.0 is
c       returned for latitude, longitude, and altitude
c
c  Inputs:
c       locatn:   character*3 location code mentioned above
c
c  Outputs:
c       lat:      r*8 latitude corresponding to locatn in decimal degrees 
c       long:     r*8 longitude corresponding to locatn in decimal degrees
c       alt:      r*8 altitude corresponding to locatn in decimal kilometers
c==================================================================
c
      subroutine xloc ( locatn, lat, long, alt )
c
      implicit none

      character locatn*(*)
c
      real*8 lat,long,alt
c
      if(index(locatn,'MCM').ne.0) then
        lat=-77.847d0
        long=166.728d0
        alt=0.1d0
      elseif (index(locatn,'JPL').ne.0) then
        lat=34.200d0
        long=-118.172d0
        alt=0.34d0
      elseif (index(locatn,'MES').ne.0) then
        lat=34.205d0
        long=-118.171d0
        alt=0.46d0
      elseif (index(locatn,'TMF').ne.0) then
        lat=34.3820d0                     ! lat=34.3816
        long=-117.6773d0                  ! long=-117.68
        alt=2.258d0                       ! alt=2.29
      elseif (index(locatn,'ARC').ne.0) then
        lat=37.43d0
        long=-122.08d0
        alt=0.01d0
      elseif (index(locatn,'PAL').ne.0) then
        lat=31.78d0
        long=-95.70d0
        alt=0.1d0
      elseif (index(locatn,'FTS').ne.0) then
        lat=34.48d0
        long=-104.22d0
        alt=1.26d0
      elseif (index(locatn,'DAG').ne.0) then
        lat=34.856d0
        long=-116.790d0
        alt=0.626d0
      elseif (index(locatn,'MTB').ne.0) then
        lat=37.584d0
        long=-118.235d0
        alt=3.80d0
      elseif (index(locatn,'LYL').ne.0) then
        lat=56.858d0
        long=-101.066d0
        alt=0.354d0
      elseif (index(locatn,'FAI').ne.0) then
        lat=64.830d0             ! 64:49:49.2  (mk4 50 m northwest of gps not
        long=-147.614d0          ! 147:36:50.4  accounted for in lat. and lon.)
        alt=0.182d0              ! 597 ft
      elseif (index(locatn,'ESN').ne.0) then
        lat=67.889d0
        long=21.085d0
        alt=0.271d0
      else
        lat=0.0d0
        long=0.0d0
        alt=0.0d0
        write(*,*)' Unknown location...',locatn
      endif
      return
      end
