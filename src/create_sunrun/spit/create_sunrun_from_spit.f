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
      integer*4 iy,im,id,hh,mm,ss,pkl,prl,
     & fnbc,lnbc,fbc,ispe,iend,dtype,nsp,nss,nip,
     & lr,lrt,ls,lunr,luns,lunt,possp,platform,istat,one,object
      parameter (lunr=14,luns=15,lunt=16,one=1)
c
      real*8 tins,pins,hins,tout,pout,hout,
     & obalt,wavtkr,oblat,oblon,
     & fsf,tcorr,nus,nue,lwn,sia,sis,aipl,
     & fxv,lxv,apt,dur,vel,phr,res
c
      character apf*2,dl*1,ext*3,spfmt*2,logfile*20,
     & outfile*64,path*128,root*64,dplist*80,specname*32,
     & string*214,user*8,col1*1
c
      write(6,*)'SUNRUN     Version 1.1.0    26-Nov-2005    GCT'
      col1=' '
c
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c     Platform specification:        DG000909
c     Also need to select path for m4head.inc in gethead.f
      call getenv('LOGNAME',user)
      if(user.ne.'        ')then
         platform=0               !0=Sun, 1=PC-Linux, 2=PC-Win32
         root='/home/toon/ggg/'
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
     &'    Nus    Nue      FSF       lasf     wavtkr   sia   sis  aipl'
c
c
      object=2 ! Sun
      oblat=78.920
      oblon=11.920
      obalt=0.100
      tcorr=0.0d0       ! Times are in UT
      fsf=0.99999999d0
      wavtkr=9900.0d0
      aipl=0.002        ! Airmass-Independent Path Length (km)
c
      do ispe=1,999999  !---------Main loop over spectra----------
         read(luns,'(a)',end=99) string 
        write(*,*)string(160:)
          read(string(160:),*)tins,pins,hins,tout,pout,hout

c  Remove any leading spaces/blanks.
         string=string(fnbc(string):)
         ls=fbc(string)
         specname=string(:ls-1)
         if(lgt(specname(1:6),'050000')) then
            nus=4000.0d0
            nue=11000.0d0
         else
            nus=6000.0d0
            nue=11000.0d0
         endif
c
c  find the spectral file, return the PATH to the spectrum
         dplist=root(:lrt)//'config'//dl//'data_part.lst'
         call gindfile(dplist,specname,path)
         if(lnbc(path).le.0) then
            write(6,*) ' Not Found : ',specname
            stop
         endif

         call read_opus_header(path,iend,dtype,nsp,fxv,lxv,iy,im,
     &   id,hh,mm,ss,apt,nss,dur,vel,apf,phr,res,lwn,nip,pkl,prl,
     &   possp,oblat,oblon,obalt,tins,pins,hins,tout,pout,hout,
     &   sia,sis)

55       call write_sunrun(lunt,col1,specname,object,tcorr,oblat,
     &   oblon,obalt,tins,pins,hins,tout,pout,hout,nus,nue,fsf,
     &   lwn,wavtkr,sia,sis,aipl,istat)
c
      end do ! -------------Main loop over spectra----------------
c
 99   close(luns)
      close(lunt)
      write(*,*) ispe-1, ' spectra found'
      stop
      end
