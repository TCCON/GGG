c  Program to create a GGG-compatible sunrun from an ascii list
c  of OPUS spectra.
c
c  To create the input file for this program run list_maker.
c
c  The sunrun program can then be immediately run from there.
c  It will search the local directory for the .gnd file, and
c  write the .gop file to the sunruns/ directory.
c
      implicit none
      include "../../ggg_int_params.f"
      include "../params.f"

      integer*4 
     & fnbc,lnbc,fbc,ispe,
     & lr,lrt,ls,
     & istat,object
c
      real*8 tins,pins,hins,tout,pout,hout,
     & obalt,wavtkr,oblat,oblon,
     & fsf,tcorr,nus,nue,lwn,sia,sis,aipl,tel_mag
c
      character 
     & spfmt*2,
     & col1*1,                    !first column of runlog record
     & dl*1,                      !forward or backward slash
     & ext*3,                     !geometry ['air','bal','gnd','lab',orb','syn']
     & gggdir*(mpath),            !ggg directory path (GGGPATH?)
     & specname*(nchar)           !spectrum name

c
      write(6,*)'create_sunrun   Version 1.1.1   9-Mar-2006   GCT'
      col1=' '
c
      call get_ggg_environment(gggdir,dl)
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
      outfile=gggdir(:lrt)//'sunruns'//dl//ext//dl//logfile(:lr-2)//'br'
      open(lunt,file=outfile,status='unknown')
      write(lunt,'(a)') ' Spectrum_File_Name   Obj   tcorr    oblat'//
     &'   oblon   obalt   tins  pins   hins   tout  pout   hout'//
     &'    Nus    Nue      FSF       lasf     wavtkr   sia   sis'
c
      object=3 ! Sun
      tcorr=0
      fsf=0.99999999d0
      fsf=1.0d0
      wavtkr=9900.0d0
      oblat=34.380d0
      oblon=-117.677
      obalt=2.270d0
      tins=25.0d0
      pins=775.8d0
      hins=20.0d0
      tout=15.0d0
      pout=775.8d0
      hout=20.0d0
      lwn=15798.0d0
      wavtkr=9900.d0
      sia=999.9d0
      sis=0.0d0
      aipl=0.003       ! Airmass-Independent Path Length (km)
      tel_mag=1.0      ! Telescope magnification (dimensionless)
c
      do ispe=1,999999  !---------Main loop over spectra----------
         read(luns,'(a)',end=99) string 

c  Remove any leading spaces/blanks.
         string=string(fnbc(string):)
         ls=fbc(string)
         specname=string(:ls-1)
         nus=0.0d0
         nue=50000.0d0
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
