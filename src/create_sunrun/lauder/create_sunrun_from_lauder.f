c  Program to create a GGG-compatible sunrun from an ascii list
c  of OPUS spectra.
c
      implicit none
      include "../../gfit/ggg_int_params.f"
      include "../params.f"

      integer*4 
     & fnbc,lnbc,fbc,ispe,lr,lrt,ls,
     & istat,object,idum
c
      real*8 tins,pins,hins,tout,pout,hout,
     & obalt,wavtkr,oblat,oblon,
     & fsf,tcorr,nus,nue,lwn,sia,sis,aipl, tel_mag
c
      character 
     & spfmt*2,
     & col1*1,                    !first column of runlog record
     & dl*1,                      !forward or backward slash
     & ext*3,                     !geometry ['air','bal','gnd','lab',orb','syn']
     & gggdir*(mpath),            !ggg directory path (GGGPATH?)
     & specname*(nchar)           !spectrum name

      idum=mfilepath ! Avoid compiler warning (unused parameter)
      idum=mauxcol  ! Avoid compiler warning (unused parameter)
      idum=mcolvav  ! Avoid compiler warning (unused parameter)
      idum=mgas     ! Avoid compiler warning (unused parameter)
      idum=mlev     ! Avoid compiler warning (unused parameter)
      idum=mrow_qc  ! Avoid compiler warning (unused parameter)
      idum=mspeci   ! Avoid compiler warning (unused parameter)
      idum=mvmode   ! Avoid compiler warning (unused parameter)
      idum=ncell    ! Avoid compiler warning (unused parameter)
      idum=nchar    ! Avoid compiler warning (unused parameter)
c
      write(6,*)'SUNRUN     Version 1.1.0    31-Oct-2005    GCT'
      col1=' '
c
      call get_ggg_environment(gggdir, dl)
      lrt=lnbc(gggdir)       !Length of root
      lr=0
      do while(lr.eq.0)
         if (iargc() == 0) then
            write(6,'(a)') 'Enter name of input file (e.g. pa2004.gnd):'
            read(*,'(a)') logfile
         elseif (iargc() == 1) then
            call getarg(1, logfile)
         else
            stop 'Usage: $gggpath/bin/create_sunrun pa2004.gnd'
         endif

         lr=lnbc(logfile)
      end do
      ext(1:1)=logfile(lr-2:lr-2)
      if(ext(1:1).eq.'a') ext(1:3)='air'
      if(ext(1:1).eq.'b') ext(1:3)='bal'
      if(ext(1:1).eq.'g') ext(1:3)='gnd'
      if(ext(1:1).eq.'l') ext(1:3)='lab'
      if(ext(1:1).eq.'o') ext(1:3)='orb'
      if(ext(1:1).eq.'s') ext(1:3)='syn'

      spfmt=logfile(lr-1:lr)

      open(luns,file=logfile,status='old')
      outfile=gggdir(:lrt)//'sunruns'//dl//ext//dl//logfile(:lr-2)//'op'
      open(lunt,file=outfile,status='unknown')
      write(lunt,'(a)') ' Spectrum_File_Name   Obj   tcorr    oblat'//
     &'   oblon   obalt   tins  pins   hins   tout  pout   hout'//
     &'    Nus    Nue      FSF       lasf     wavtkr   sia   sis'//
     &'   AIPL   TM'
c
      object=2 ! Sun
c set tcorr once parsed filename (as file-dependent) 
c      tcorr=-(12*3600.0d0) ! 12 hours ahead of UT
      oblat=-45.050d0
      oblon=169.680d0
      obalt=0.370d0
      tins=25.0d0
      pins=1.0d0     ! Revised in an external append of met data 
      hins=20.0d0
      tout=0.0d0     ! Replaced in an external append of met data 
      pout=0.0d0     ! Replaced in an external append of met data 
      hout=50.0d0
c      fsf=0.9999999d0 ! set instrument dependent fsf below
      lwn=15798.0d0
      wavtkr=9900.d0
      sia=999.9d0
      sis=0.0d0
      aipl=0.002       ! Airmass-Independent Path Length (km)
      tel_mag=1.0      
c
      do ispe=1,999999  !---------Main loop over spectra----------
         read(luns,'(a)',end=99) string 

c  Remove any leading spaces/blanks.
         string=string(fnbc(string):)
         ls=fbc(string)
         specname=string(:ls-1)
         if (index(specname,'a.') .gt. 0) then
             !times for spectra transformed using opus-ipp are in UTC
             tcorr=(0*3600.0d0) ! UTC for spectra transformed using opus-ipp
         else
             !times for Bruker interferograms and spectra are UTC+12
             tcorr=-(12*3600.0d0) ! 12 hours ahead of UT
         endif
         if (specname(1:1) .eq. 'l') fsf=1.000000d0 ! 120HR
         if (specname(1:1) .eq. 'L') fsf=0.999995d0 ! 120M
         if(specname(2:2).eq.'w' .or. specname(2:2).eq.'W') then   ! InGaAs det.
            nus=4000.0d0
            nue=8500.0d0
         elseif(specname(2:2).eq.'h') then   ! InGaAs det. new filename convention
            nus=4000.0d0
            nue=8500.0d0
         elseif(specname(2:2).eq.'v' .or. specname(2:2).eq.'V') then   ! Si det.
            nus=12950.0d0
            nue=14800.0d0
         elseif(specname(16:16).eq.'c') then   ! InSb detector
            nus=1800.0d0
            nue=6000.0d0
         else
            write(*,*) 'Unknown spectrum type'
            stop
         endif
c
c  find the spectral file, return the PATH to the spectrum
          dplist=gggdir(:lrt)//'config'//dl//'data_part.lst'
          call gindfile(dplist,specname,path)
          if(lnbc(path).le.0) then
            write(6,*) ' Not Found: ',specname
            stop
          endif

          call write_sunrun(lunt,col1,specname,object,tcorr,oblat,
     &    oblon,obalt,tins,pins,hins,tout,pout,hout,nus,nue,fsf,
     &    lwn,wavtkr,sia,sis,aipl,tel_mag,istat)
c
      end do ! -------------Main loop over spectra----------------
c
 99   close(luns)
      close(lunt)
      write(*,*) ispe-1, ' spectra found'
      stop
      end

