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
     & fnbc,lnbc,fbc,ispe,iend,dtype,nsp,nip,
     & lr,lrt,ls,possp,istat,object
c
      real*8 tins,pins,hins,tout,pout,hout,
     & wavtkr,oblat,oblon,obalt,lfl,hfl,foc,dfr,
     & fsf,tcorr,nus,nue,lwn,sia,sis,aipl,tel_mag,
     & fxv,lxv,apt,dur,vel,phr,res
c
      character 
     & col1*1,                    !first column of runlog record
     & apf*2,                     !apodization function (e.g. BX N2, etc)
     & dl*1,                      !forward or backward slash
     & ext*3,                     !geometry ['air','bal','gnd','lab',orb','syn']
     & gggdir*(mpath),            !ggg directory path (GGGPATH?)
     & specname*(nchar)           !spectrum name

c
      write(6,*)'SUNRUN     Version 1.1.1    27-Nov-2007    GCT'
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
c      write(*,*)logfile
      ext(1:1)=logfile(lr-2:lr-2)
      if(ext(1:1).eq.'a') ext(1:3)='air'
      if(ext(1:1).eq.'b') ext(1:3)='bal'
      if(ext(1:1).eq.'g') ext(1:3)='gnd'
      if(ext(1:1).eq.'l') ext(1:3)='lab'
      if(ext(1:1).eq.'o') ext(1:3)='orb'
      if(ext(1:1).eq.'s') ext(1:3)='syn'

      open(luns,file=logfile,status='old')
      outfile=gggdir(:lrt)//'sunruns'//dl//ext//dl//logfile(:lr-2)//'op'
      open(lunt,file=outfile,status='unknown')
      write(lunt,'(a)') ' Spectrum_File_Name   Obj   tcorr    oblat'//
     &'   oblon   obalt   tins  pins   hins   tout  pout   hout'//
     &'    Nus    Nue      FSF       lasf     wavtkr   sia   sis  aipl'
c
c
      object=2 ! Sun
      oblat=53.
      oblon=7
      obalt=0.027
      tcorr=0.0d0       ! Times are in UT
      fsf=1.0d0
      wavtkr=9900.0d0
      aipl=0.004        ! Airmass-Independent Path Length (km)
c
      do ispe=1,999999  !---------Main loop over spectra----------
         read(luns,'(a)',end=99) string 
c        write(*,*)string(160:)
c          read(string(160:),*)tins,pins,hins,tout,pout,hout

c  Remove any leading spaces/blanks.
         string=string(fnbc(string):)
         ls=fbc(string)
         specname=string(:ls-1)
         if(lgt(specname(1:6),'050000')) then
            nus=4000.0d0
            nue=15000.0d0
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

c         write(*,*)'Calling read_opus_header:',path
         call read_opus_header(path,iend,dtype,nsp,fxv,lxv,iy,im,
     & id,hh,mm,ss,ms,apt,dur,vel,apf,phr,res,lwn,foc,nip,dfr,
     & pkl,prl,gfw,gbw,lfl,hfl,possp,oblat,oblon,obalt,tins,pins,hins,
     & tout,pout,hout,sia,sis)
c         write(*,*)'Called read_opus_header'

         call write_sunrun(lunt,col1,specname,object,tcorr,oblat,
     &   oblon,obalt,tins,pins,hins,tout,pout,hout,nus,nue,fsf,
     &   lwn,wavtkr,sia,sis,aipl,tel_mag,istat)
c
      end do ! -------------Main loop over spectra----------------
c
 99   close(luns)
      close(lunt)
      write(*,*) ispe-1, ' spectra found'
      stop
      end
