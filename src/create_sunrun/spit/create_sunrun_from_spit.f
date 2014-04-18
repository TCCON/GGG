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

      integer*4 iy,im,id,hh,mm,ss,ms,pkl,prl,gfw,gbw,
     & fnbc,lnbc,fbc,ispe,iend,dtype,nsp,nip,dfr,
     & lr,lrt,ls,
     & possp,istat,object,lst
c
      real*8 tins,pins,hins,tout,pout,hout,
     & wavtkr,oblat,oblon,obalt,lfl,hfl,foc,
     & fsf,tcorr,nus,nue,lwn,sia,sis,vdc,aipl,
     & fxv,lxv,apt,dur,vel,phr,res,lse,lsu
c
      character
     & spfmt*2,
     & col1*1,                    !first column of runlog record
     & apf*2,                     !apodization function (e.g. BX N2, etc)
     & dl*1,                      !forward or backward slash
     & ext*3,                     !geometry ['air','bal','gnd','lab',orb','syn']
     & gggdir*(mpath),            !ggg directory path (GGGPATH?)
     & specname*(nchar)           !spectrum name

c
      write(6,*)'SUNRUN     Version 1.1.1    28-Nov-2007    GCT'
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
     &'    Nus    Nue      FSF       lasf     wavtkr   sia   sis  aipl'
c
c
      object=2 ! Sun
      oblat=78.920
      oblon=11.920
      obalt=0.010
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
         dplist=gggdir(:lrt)//'config'//dl//'data_part.lst'
         call gindfile(dplist,specname,path)
         if(lnbc(path).le.0) then
            write(6,*) ' Not Found : ',specname
            stop
         endif

       call read_opus_header(path,iend,dtype,nsp,fxv,lxv,iy,im,
     &  id,hh,mm,ss,ms,apt,dur,vel,apf,phr,res,lwn,foc,nip,dfr,
     &  pkl,prl,gfw,gbw,lfl,hfl,possp,oblat,oblon,obalt,
     &  tins,pins,hins,tout,pout,hout,sia,sis,vdc,lst,lse,lsu)

         call write_sunrun(lunt,col1,specname,object,tcorr,oblat,
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
