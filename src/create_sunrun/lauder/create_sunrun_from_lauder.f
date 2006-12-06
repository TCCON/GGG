c  Program to create a GGG-compatible sunrun from an ascii list
c  of OPUS spectra.
c
c  To create the input file for this program go to the
c  directory containing the relevant spectra and create a list
c       ls -1 db20051024* > temp_file
c  Unfortunately this will put all the InGaAs spectra in the
c  first half of the file, and all the Si in the second half.
c
c  So sort it so that it is in chronological order.
c       sort -k 1.3,1.10 -k 1.18 temp_file > db20051018.gnd
c
c  The .gnd file is intended to reside in the local directory,
c  which is usually the directory containing the spectra.
c  The sunrun program can then be immediately run from there.
c  It will search the local directory for the .gnd file, and
c  write the .gop file to the sunruns/ directory.
c
      integer*4 
     & fnbc,lnbc,fbc,ispe,
     & lr,lrt,ls,lunr,luns,lunt,platform,istat,one,object
      parameter (lunr=14,luns=15,lunt=16,one=1)
c
      real*8 tins,pins,hins,tout,pout,hout,
     & obalt,wavtkr,oblat,oblon,
     & fsf,tcorr,nus,nue,lwn,sia,sis,aipl,tel_mag
c
      character dl*1,ext*3,spfmt*2,logfile*20,
     & outfile*64,path*128,root*64,dplist*80,specname*32,
     & string*163,user*8,col1*1
c
      write(6,*)'create_sunrun   Version 1.1.1   9-Mar-2006   GCT'
      col1=' '
c
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c     Platform specification:        DG000909
c     Also need to select path for m4head.inc in gethead.f
      call getenv('LOGNAME',user)
      if(user.ne.'        ')then
         platform=0               !0=Sun, 1=PC-Linux, 2=PC-Win32
c         root='/home/toon/ggg/'
         call getenv('GGGPATH',root)
         dl='/'
         root=root(:lnbc(root))//dl
      else
         platform=2               !0=Sun, 1=PC-Linux, 2=PC-Win32
         dl=char(92)  ! backslash ('\')
         root='g:'//dl
         user='PC-Win'
      endif
      lrt=lnbc(root)       !Length of root
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      lr=0
      do while(lr.eq.0)
         write(6,'(a,$)') 'Enter name of input file (e.g. pa2004.gnd): '
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
      outfile=root(:lrt)//'sunruns'//dl//ext//dl//logfile(:lr-2)//'op'
      open(lunt,file=outfile,status='unknown')
      write(lunt,'(a)') ' Spectrum_File_Name   Obj   tcorr    oblat'//
     &'   oblon   obalt   tins  pins   hins   tout  pout   hout'//
     &'    Nus    Nue      FSF       lasf     wavtkr   sia   sis'
c
      object=2 ! Sun
      tcorr=-(12*3600.0d0) ! 12 hours ahead of UT
      fsf=0.99999999d0
      wavtkr=9900.0d0
      oblat=-45.050d0
      oblon=169.680d0
      obalt=0.370d0
      tins=25.0d0
      pins=1.0d0     ! 120 HR
      hins=20.0d0
      tout=15.0d0
      pout=970.0d0
      hout=50.0d0
      fsf=0.9999999d0
      lwn=15798.0d0
      wavtkr=9900.d0
      sia=999.9d0
      sis=0.0d0
      aipl=0.002       ! Airmass-Independent Path Length (km)
      tel_mag=1.0      ! Telescope magnification (dimensionless)
c
      do ispe=1,999999  !---------Main loop over spectra----------
         read(luns,'(a)',end=99) string 

c  Remove any leading spaces/blanks.
         string=string(fnbc(string):)
         ls=fbc(string)
         specname=string(:ls-1)
         if(specname(16:16).eq.'a') then      ! InGaAs detector
            nus=4000.0d0
            nue=11000.0d0
         elseif(specname(16:16).eq.'b') then  ! Si detector
            if(specname(15:15).eq.'0') then
               nus=11000.0d0  ! No Red Filter
               nue=13500.0d0  ! No Red Filter
            else
               nus=11000.0d0  ! Red Filter
               nue=15000.0d0  ! Red Filter
            endif
         elseif(specname(16:16).eq.'c') then   ! InSb detector
            nus=1800.0d0
            nue=6000.0d0
         else
            write(*,*) 'Unknown spectrum type'
            stop
         endif
c
c  find the spectral file, return the PATH to the spectrum
          dplist=root(:lrt)//'config'//dl//'data_part.lst'
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
