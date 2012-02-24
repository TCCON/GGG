c  Program to create a GGG-compatible sunrun from an ascii list
c  of OPUS spectra.
c
      implicit none
      include "../../ggg_int_params.f"
      include "../params.f"

      integer*4 
     & fnbc,lnbc,fbc,ispe,lr,lrt,ls,
     & istat,object
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

c
      write(6,*)'create_sunrun     Version 1.1.1    24-Apr-2008    GCT'
      col1=' '
c
      call get_ggg_environment(gggdir, dl)
      lrt=lnbc(gggdir)       !Length of root
      lr=0
      do while(lr.eq.0)
         write(6,'(a)') 'Enter name of input file (e.g. pa2004.gnd): '
         read(*,'(a)') logfile
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
      tcorr=(0*3600.0d0) ! UT
      oblat=28.3079d0
      oblon=-16.5000d0
      obalt=2.370d0
      tins=25.0d0
c      pins=770.0d0     ! Revised in an external append of met data 
      hins=20.0d0
c      tout=0.0d0     ! Replaced in an external append of met data 
c      pout=0.0d0     ! Replaced in an external append of met data 
c      hout=50.0d0
      fsf=0.99999999d0 ! set instrument dependent fsf below
      lwn=15798.0d0
      wavtkr=9900.d0
      sia=999.9d0
      sis=0.0d0
      aipl=0.002       ! Airmass-Independent Path Length (km)
      tel_mag=1.0      
      nus=3800.0d0
      nue=8500.0d0
c
      do ispe=1,999999  !---------Main loop over spectra----------
         read(luns,'(a)',end=99) string 

c  Remove any leading spaces/blanks.
         string=string(fnbc(string):)
         ls=fbc(string)
         specname=string(:ls-1)

         if(specname(1:6).eq.'070520') then
           tout=9.3
           pout=764.7
           hout=40.0
         elseif(specname(1:6).eq.'070830') then
           tout=19.4
           pout=772.6
           hout=35.0
         elseif(specname(1:6).eq.'071202') then
           tout=6.5
           pout=771.1
           hout=45.0
         endif
         pins=pout
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

