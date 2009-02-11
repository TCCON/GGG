c  Program to create a GGG-compatible runlog from the sunrun file.
c  which is created by running create_sunrun_from_xxxx.
c
      implicit none

      integer*4 bytepw,ifirst,ilast,iy,im,id,jj,jd,i,
     & lnbc,ispe,ifmin,ifmax,istat,nlhead,ncol,
     & lr,lrt,luns,lunt,possp,doy,one,object
      parameter (luns=15,lunt=16,one=1)
c
      real*8 amal,fovi,fovo,gmt,tins,pins,hins,tout,pout,hout,
     & obalt,asza,azim,snr,wavtkr,zoff,zpoff,oblat,oblon,opd,
     & lasf,fmin,fmax,fsf,delwav,tcorr,fvsi,wspd,wdir,aipl,
     & tel_mag,eorv,ervc,osds,tplat,tplon,tpalt,site_solar_noon
c
c      logical flexst
      character apf*2,dl*1,ext*3,spfmt*2,logfile*40,outfile*64,
     & path*128,root*64,dplist*80,specname*35,col1*1,version*52
c
      version=' CREATE_RUNLOG    Version 8.5.1    05-Feb-2009   GCT'
      col1=' '
      iy=0
      im=0
      id=30
      amal=0.0d0
c
ccxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cc   Root path specification
c      call getenv('LOGNAME',user)
c      if(user.ne.'        ')then     !Sun, PC-Linux
c         call getenv('GGGPATH',root)
c         dl='/'
c         root=root(:lnbc(root))//dl
c      else                           !PC-Win32
c         dl=char(92)  ! backslash ('\')
c         root='g:'//dl
c         user='PC-Win'
c      endif
c      lrt=lnbc(root)       !Length of root
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      dl='/'
      call getenv('GGGPATH',root)
      root=root(:lnbc(root))//dl
      lrt=lnbc(root)       !Length of root

      lr=0
      do while(lr.eq.0)
         write(6,'(a)') 'Enter sunrun (e.g. fts89avg.gop): '
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
c      write(*,*)root(:lrt)//'sunruns'//dl//ext//dl//logfile
      open(luns,file=root(:lrt)//'sunruns'//dl//ext//dl//logfile,
     $status='old')
      read(luns,*) nlhead, ncol
      do i=2,nlhead
         read(luns,*)
      end do
      outfile=root(:lrt)//'runlogs'//dl//ext//dl//logfile(:lr-2)//'rl'
c
      open(lunt,file=outfile,status='unknown')
      write(lunt,*)3,35
      write(lunt,'(a)') version
      write(lunt,'(a)')
     &  '  Spectrum_File_Name                 Year  Day  Hour'//
     &  '   oblat    oblon   obalt    ASZA   POFF    AZIM   OSDS'//
     &  '    OPD   FOVI  FOVO'//
     &  '  AMAL  IFIRST   ILAST    DELTA_NU   POINTER  BPW ZOFF SNR'//
     &  '  APF tins  pins  hins   tout   pout  hout'//
     &  '  fvsi  wspd  wdir  lasf    wavtkr  aipl'
c
      do ispe=1,9999999  !---------Main loop over spectra----------
         call read_sunrun(luns,col1,specname,object,tcorr,oblat,
     &   oblon,obalt,tins,pins,hins,tout,pout,hout,fvsi,wspd,wdir,
     &   fmin,fmax,fsf,lasf,wavtkr,aipl,tel_mag,istat)
         if(istat.ne.0) go to 99
c
c  find the spectral file, return the PATH to the spectrum
         dplist=root(:lrt)//'config'//dl//'data_part.lst'
         call gindfile(dplist,specname,path)
         if(lnbc(path).le.0) then
            write(6,*) ' Not Found: "'//specname//'"'
            stop
         endif

c         write(*,*)'create_runlog: ',specname,path
         call rdsphead(spfmt,specname,path,ifirst,ilast,possp,
     &   bytepw,apf,delwav,opd,fovi,snr,
     &   iy,im,id,gmt,lasf,wavtkr)

         fovo=fovi/tel_mag
c         write(*,*) path,ifirst,ilast,possp,bytepw,apf,delwav

c  Fudge the starting index to prevent use of noise regions.
c  Multiply by delwav to allow for reverse-ordered spectra
c  in which IFIRST, ILAST & DELWAV are all negative.
         ifmin=nint(fmin/delwav)
         if(ifmin*delwav.gt.ifirst*delwav) then
            possp=possp+iabs(bytepw)*(ifmin-ifirst)
            ifirst=ifmin
         endif
c
c  Fudge the ending index to prevent use of noise regions.
         ifmax=nint(fmax/delwav)
c         write(*,*)ifmax,ilast,fmax,ilast*delwav
         if(ifmax*delwav.lt.ilast*delwav) then
            ilast=ifmax
         endif
c
         delwav=delwav*fsf
         gmt=gmt+tcorr/3600.0d0
c
c  Calculate Day of Year (DOY)
         call julian(iy,im,id,jd)
         call julian(iy,one,one,jj)
         doy=jd-jj+1
c
c  We want all observations to have the day number of local noon,
c  even though the UT day number may be different. The following
c  makes all spectra acquired on the same local day have the same
c  day number. To achieve this, the UT time must be allowed to go
c  -ve or to exceed 24.
         site_solar_noon=12-oblon/15  !  ut time of local noon
         if((gmt-site_solar_noon).gt.+12) then
             gmt=gmt-24.0
             doy=doy+1
         elseif((gmt-site_solar_noon).lt.-12) then
             gmt=gmt+24.0
             doy=doy-1
         endif
c         write(*,*)path(:lnbc(path))
c
         call zenaz(object,oblat,oblon,obalt,iy,one,doy,
     &   gmt/24.d0,asza,azim,eorv,ervc,tplat,tplon,tpalt)
         osds=1.e+06*(eorv+ervc)/3.d+8  ! Observer-Sun Doppler Stretch (ppm)

         zoff=0.0
         zpoff=0.0
         call write_runlog(lunt,col1,specname,iy,doy,gmt,oblat,
     &   oblon,obalt,asza,zpoff,azim,osds,opd,fovi,fovo,amal,
     &   ifirst,ilast,
     &   delwav,possp,bytepw,zoff,snr,apf,tins,pins,hins,tout,
     &   pout,hout,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)
         if(mod(ispe,1000).eq.0) write(*,*) ispe
c
      end do ! -------------Main loop over spectra----------------
c
 99   close(luns)
      close(lunt)
      write(*,*) ispe-1, ' spectra found'
      stop
      end

