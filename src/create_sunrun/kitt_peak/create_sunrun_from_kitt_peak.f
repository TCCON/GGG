c  Program to create a GGG-compatible sunrun from an ascii list of spectra.
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

      integer*4 iy,im,id,hh,mm,ss,
     & fnbc,lnbc,fbc,ispe,iend,dtype,
     & lr,lrt,ls,
     & possp,istat,object
c
      real*8 tins,pins,hins,tout,pout,hout,
     & obalt,wavtkr,oblat,oblon,aipl,tel_mag,
     & fsf,tcorr,nus,nue,sia,fvsi,wspd,wdir,lwn

      integer*4 nsp,nip,pkl,prl

      real*8 vel, phr, res, apt,  dur, lasf

      character 

     & col1*1,                    !first column of runlog record
     & apf*2,                     !apodization function (e.g. BX N2, etc)
     & dl*1,                      !forward or backward slash
     & ext*3,                     !geometry ['air','bal','gnd','lab',orb','syn']
     & gggdir*(mpath),            !ggg directory path (GGGPATH?)
     & specname*(nchar),          !spectrum name
     & version*64                 !current program version
c
      version=
     &' create_sunrun_from_kp     Version 1.2.2    18-Dec-2011    GCT'
      write(6,*) version
c
      col1=' '
      call get_ggg_environment(gggdir, dl)
      lrt=lnbc(gggdir)               ! Length of root
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

      open(luns,file=logfile,status='old')
      outfile=gggdir(:lrt)//'sunruns'//dl//ext//dl//logfile(:lr-2)//'fi'
      open(lunt,file=outfile,status='unknown')
      write(lunt,*)3,23
      write(lunt,'(a)') version
      write(lunt,'(a)')
     &' Spectrum_File_Name                                      '//
     &'  Obj   tcorr'//
     &'   oblat   oblon   obalt   tins   pins   hins'//
     &'   tout  pout   hout   sia   fvsi  wspd   wdir'//
     &'   Nus    Nue      FSF      lasf     wavtkr'//
     &'   AIPL   TM'
c
      object=2          ! Sun=2; Moon=1
      tcorr=7*3600.0d0  ! Kitt peak is 7 hours behind UT
      wavtkr=9900.0d0
      fsf=1.00000000d0
      aipl=0.226        ! Airmass-Independent Path Length (km)
      tel_mag=50.       ! Telescope Magnification

      do ispe=1,999999  !---------Main loop over spectra----------
         read(luns,'(a)',end=99) string 

c  Remove any leading spaces/blanks.
         string=string(fnbc(string):)
         ls=fbc(string)
         specname=string(:ls-1)
c
c  find the spectral file, return the PATH to the spectrum
         dplist=gggdir(:lrt)//'config'//dl//'data_part.lst'
         call gindfile(dplist,specname,path)
         if(lnbc(path).le.0) then
            write(6,*) ' Not Found: ',specname
            stop
         endif

      write(*,*)'Calling read_fits_header...',path
      call read_fits_header(path,iend,dtype,nsp,nus,nue,iy,im,id,
     & hh,mm,ss,apt,dur,vel,apf,phr,res,lasf,nip,pkl,prl,
     & possp,oblat,oblon,obalt,tins,pins,hins,tout,pout,hout)

      write(*,*)'Called read_fits_header.'
c YZH frequency shift adjustments:
      fsf=1.0d0
      if(iy.lt.1978) then
         fsf=0.999999d0
      else if(iy.eq.1978.and.im.eq.5) then
         fsf=0.999999d0
      else if(iy.eq.1978.and.im.eq.9) then
         fsf=0.999998d0
      else if(iy.eq.1978.and.im.eq.11) then
         fsf=0.9999965d0
         fsf=0.9999985d0
      else if(iy.eq.1979.and.im.le.8) then
         fsf=0.999999d0
      else if(iy.eq.1979.and.im.eq.10) then
         fsf=0.999997d0
      else if(iy.eq.1979.and.(im.eq.11.or.im.eq.12)) then
         fsf=0.999999d0
      else if(iy.eq.1980.and.im.eq.2) then
         fsf=0.999997d0
      else if(iy.eq.1980.and.im.eq.3) then
         fsf=0.999996d0
      else if(iy.eq.1980.and.im.eq.5) then
         fsf=0.999998d0
      else if(iy.eq.1980.and.im.eq.6) then
         fsf=0.999999d0
      else if(iy.eq.1980.and.im.eq.7) then
         fsf=0.999997d0
      else if(iy.eq.1980.and.im.eq.8) then
         fsf=0.999997d0
      else if(iy.eq.1980.and.im.eq.9) then
         fsf=0.999999d0
c      else if(iy.eq.1980.and.(im.le.3.or.im.eq.7.or.im.eq.8)) then
c         fsf=0.999996d0
      else if(iy.eq.1981.and.im.eq.1) then
         fsf=0.999999d0
      else if(iy.eq.1981.and.im.eq.2) then
         fsf=0.999995d0
      else if(iy.eq.1981.and.im.eq.5) then
         fsf=1.000001d0
      else if(iy.eq.1981.and.im.eq.7) then
         fsf=0.999998d0
      else if(iy.eq.1981.and.im.eq.8) then
         fsf=0.999999d0
      else if(iy.eq.1981.and.im.eq.9) then
         fsf=0.999999d0
      else if(iy.eq.1981.and.im.eq.10) then
         fsf=0.999998d0
      else if(iy.eq.1981.and.im.eq.11) then
         fsf=0.999999d0
      else if(iy.eq.1982.and.im.eq.10.and.id.eq.1) then
         fsf=1.000001d0
      else if(iy.eq.1982.and.(im.ge.3.and.im.le.11)) then
         fsf=0.999999d0
      else if(iy.eq.1983.and.im.le.2) then
         fsf=0.9999995d0
      else if(iy.eq.1983.and.im.eq.4) then
         fsf=0.999995d0
      else if(iy.eq.1983.and.im.eq.6) then
         fsf=0.999999d0
      else if((iy.eq.1983.and.im.eq.5.and.id.eq.4)) then
         fsf=0.999995d0
      else if(iy.eq.1985.and.im.eq.12.and.id.eq.13) then
         fsf=1.000001d0
      else if(iy.eq.1990.and.im.eq.4) then
         fsf=1.000001d0
      else if(iy.eq.1991.and.im.eq.12) then
         fsf=0.999999d0
      else if(iy.eq.1995.and.im.eq.4) then
         fsf=1.000001d0
         fsf=1.0000005d0
      else if(iy.eq.1997.and.im.eq.11) then
         fsf=0.9999997d0
      endif
c
c      if(iy.gt. 1985 .or. iy .lt. 1979) then
c        fsf=0.999999d0
c      else
c        fsf=0.999996d0
c      endif

       write(*,*)'Calling write_sunrun...',specname
       fvsi=0.0d0
       sia=0.0d0
       wspd=0.0d0
       wdir=0.0d0
       call write_sunrun(lunt,col1,specname,object,tcorr,oblat,
     & oblon,obalt,tins,pins,hins,tout,pout,hout,sia,fvsi,wspd,wdir,
     & nus,nue,fsf,lwn,wavtkr,aipl,tel_mag,istat)

c       write_sunrun(lun,col1,specname,obj,tcorr,oblat,
c     &   oblon,obalt,tins,pins,hins,tout,pout,hout,sia,fvsi,
c     &   wspd,wdir,nus,nue,fsf,lasf,wavtkr,aipl,tel_mag,istat)

       write(*,*)'Called write_sunrun.'
c
      end do ! -------------Main loop over spectra----------------
c
 99   close(luns)
      close(lunt)
      write(*,*) ispe-1, ' spectra found'
      stop
      end

