c  Program to create a GGG-compatible sunrun from an ascii list of spectra.
c
      implicit none

      integer*4 iy,im,id,hh,mm,ss,
     & fnbc,lnbc,fbc,ispe,iend,dtype,
     & lr,lrt,ls,lunr,luns,lunt,possp,platform,istat,one,object
      parameter (lunr=14,luns=15,lunt=16,one=1)
c
      real*8 tins,pins,hins,tout,pout,hout,
     & obalt,wavtkr,oblat,oblon,aipl,
     & fsf,tcorr,nus,nue,sia,sis

      integer*4 nss,nsp,nip,pkl,prl

      real*8 vel, phr, res, apt,  dur, lasf

      character apf*2,dl*1,ext*3,spfmt*2,logfile*20,outfile*64,
     & path*128,root*64,dplist*80,specname*32,
     & string*163,user*8,col1*1
c
      write(6,*)'create_sunrun    Version 1.0.0    4-Nov-2005    GCT'
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

c
      spfmt=logfile(lr-1:lr)

      open(luns,file=logfile,status='old')
      outfile=root(:lrt)//'sunruns'//dl//ext//dl//logfile(:lr-2)//'fi'
      open(lunt,file=outfile,status='unknown')
      write(lunt,'(a)') '  Spectrum_File_Name   Obj   tcorr'//
     & '    oblat   oblon   obalt  tins   pins   hins '//
     & '  tout   pout  hout    Nus    Nue     FSF       lasf'//
     & '     wavtkr   sia   sis'
c
      do ispe=1,999999  !---------Main loop over spectra----------
         read(luns,'(a)',end=99) string 

c  Remove any leading spaces/blanks.
         string=string(fnbc(string):)
         ls=fbc(string)
         specname=string(:ls-1)
c
c  find the spectral file, return the PATH to the spectrum
         dplist=root(:lrt)//'config'//dl//'data_part.lst'
         call gindfile(dplist,specname,path)
         if(lnbc(path).le.0) then
            write(6,*) ' Not Found: ',specname
            stop
         endif

      call read_fits_header(path,iend,dtype,nsp,nus,nue,iy,im,id,
     & hh,mm,ss,apt,nss,dur,vel,apf,phr,res,lasf,nip,pkl,prl,
     & possp,oblat,oblon,obalt,tins,pins,hins,tout,pout,hout)

      object=2          ! Sun=2; Moon=1
      tcorr=7*3600.0d0  ! Kitt peak is 7 hours behind UT
      wavtkr=9900.0d0
      fsf=1.00000000d0
      aipl=0.226        ! Airmass-Independent Path Length (km)
c YZH frequency shift adjustments:
c      if(iy.lt.1979) then
c         fsf=0.999998d0
c      else if(iy.eq.1979.and.im.eq.5) then
c         fsf=0.9999998d0
c      else if(iy.eq.1979.and.im.le.8) then
c         fsf=0.999998d0
c      else if(iy.eq.1979.and.im.eq.10) then
c         fsf=0.999996d0
c      else if(iy.eq.1979.and.(im.eq.11.or.im.eq.12)) then
c         fsf=0.999998d0
c      else if(iy.eq.1980.and.(im.le.3.or.im.eq.7.or.im.eq.8)) then
c         fsf=0.999996d0
c      else if(iy.eq.1980.and.im.eq.6) then
c         fsf=0.999998d0
c      else if(iy.eq.1980.and.im.eq.11) then
c         fsf=0.999999d0
c      else if(iy.eq.1981.and.im.eq.1) then
c         fsf=0.999998d0
c      else if(iy.eq.1981.and.im.eq.2) then
c         fsf=0.999995d0
c      else if(iy.eq.1981.and.(im.ge.3.and.im.le.5)) then
c         fsf=0.999999d0
c      else if(iy.eq.1980.and.im.eq.5.and.id.gt.10) then
c         fsf=0.999999d0
c      else if(iy.eq.1980.and.im.eq.5.and.id.gt.10) then
c         fsf=0.999998d0
c      else if(iy.eq.1981.and.((im.ge.8.and.im.le.9).or.im.eq.11)) then
c         fsf=0.999998d0
c      else if(iy.eq.1981.and.im.eq.12) then
c         fsf=0.999999d0
c      else if(iy.eq.1982.and.(im.eq.2.or.im.eq.12)) then
c         fsf=0.999999d0
c      else if(iy.eq.1982.and.(im.ge.3.and.im.le.11)) then
c         fsf=0.999998d0
c      else if(iy.eq.1982.and.im.eq.10.and.id.eq.1) then
c         fsf=0.999999d0
c      else if(iy.eq.1983.and.im.le.2) then
c         fsf=0.9999985d0
c      else if(iy.eq.1983.and.im.eq.4) then
c         fsf=0.999993d0
c      else if((iy.eq.1983.and.im.ge.5).or.(iy.gt.1983)) then
c         fsf=0.999999d0
c      else if(iy.eq.1985.and.im.eq.12.and.id.eq.13) then
c         fsf=0.9999998d0
c      else if(iy.eq.1990.and.im.eq.4) then
c         fsf=0.9999998d0
c      else if(iy.eq.1990.and.im.eq.4) then
c         fsf=0.9999998d0
c      else if(iy.eq.1995.and.im.eq.4) then
c         fsf=0.9999997d0
c      else if(iy.eq.1997.and.im.eq.11) then
c         fsf=0.9999997d0
c      else
c         fsf=0.999997d0
c      endif
c
c      if(iy.gt. 1985 .or. iy .lt. 1979) then
c        fsf=0.999999d0
c      else
c        fsf=0.999996d0
c      endif

55    call write_sunrun(lunt,col1,specname,object,tcorr,oblat,
     & oblon,obalt,tins,pins,hins,tout,pout,hout,nus,nue,fsf,
     & lasf,wavtkr,sia,sis,aipl,istat)
c
      end do ! -------------Main loop over spectra----------------
c
 99   close(luns)
      close(lunt)
      write(*,*) ispe-1, ' spectra found'
      stop
      end

