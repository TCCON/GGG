c  Program to create a GGG-compatible runlog from the sunrun file.
c  which is created by running create_sunrun_from_xxxx.
c
      implicit none

      integer*4 bytepw,ifirst,ilast,iy,im,id,jj,jd,
     & lnbc,ispe,ifmin,ifmax,istat,
     & lr,lrt,luns,lunt,possp,doy,one,object
      parameter (luns=15,lunt=16,one=1)
c
      real*8 amal,fovi,fovo,gmt,tins,pins,hins,tout,pout,hout,
     & obalt,asza,azim,snr,wavtkr,zoff,zpoff,oblat,oblon,opd,
     & lasf,fmin,fmax,fsf,delwav,tcorr,sia,sis,aipl,tel_mag,
     & eorv,ervc,tplat,tplon,tpalt,site_solar_noon,wlimit
c
c      logical flexst
      character apf*2,dl*1,ext*3,spfmt*2,logfile*40,outfile*64,
     & path*128,root*64,dplist*80,specname*32,user*8,col1*1
c
      write(6,*)'CREATE_RUNLOG   Version 8.3.1   12-Aug-2008   GCT'
      col1=' '
      iy=0
      im=0
      id=30
      amal=0.0d0
c
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c   Root path specification
      call getenv('LOGNAME',user)
      if(user.ne.'        ')then     !Sun, PC-Linux
         call getenv('GGGPATH',root)
         dl='/'
         root=root(:lnbc(root))//dl
      else                           !PC-Win32
         dl=char(92)  ! backslash ('\')
         root='g:'//dl
         user='PC-Win'
      endif
      lrt=lnbc(root)       !Length of root
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      lr=0
      do while(lr.eq.0)
         write(6,'(a)') 'Enter sunrun (e.g. fts89avg.bal): '
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

      open(luns,file=root(:lrt)//'sunruns'//dl//ext//dl//logfile,
     $status='old')
      read(luns,*)
      outfile=root(:lrt)//'runlogs'//dl//ext//dl//logfile(:lr-2)//'rl'
c
      open(lunt,file=outfile,status='unknown')
      write(lunt,'(a)') ' Spectrum_File_Name    Year  Day  Hour'//
     &  '   oblat    oblon   obalt    ASZA   POFF    OPD   FOVI  FOVO'//
     &  '  AMAL  IFIRST   ILAST    DELTA_NU   POINTER  BPW ZOFF SNR'//
     &  '  APF tins  pins  hins   tout   pout  hout   lasf    wavtkr'//
     &  '  sia   sis   aipl'
c
      do ispe=1,9999999  !---------Main loop over spectra----------
         call read_sunrun(luns,col1,specname,object,tcorr,oblat,
     &   oblon,obalt,tins,pins,hins,tout,pout,hout,fmin,fmax,fsf,
     &   lasf,wavtkr,sia,sis,aipl,tel_mag,istat)
         if(istat.ne.0) go to 99
c
c  find the spectral file, return the PATH to the spectrum
         dplist=root(:lrt)//'config'//dl//'data_part.lst'
         call gindfile(dplist,specname,path)
         if(lnbc(path).le.0) then
            write(6,*) ' Not Found: "'//specname//'"'
            stop
         endif

         call rdsphead(spfmt,specname,path,ifirst,ilast,possp,
     &   bytepw,apf,delwav,opd,fovi,snr,
     &   iy,im,id,gmt,lasf,wavtkr)
         fovo=fovi/tel_mag
c         write(*,*) path,ifirst,ilast,possp,bytepw,apf,delwav

c  Fudge the starting index to prevent use of noise regions.
         ifmin=nint(fmin/delwav)
         if(ifmin.gt.ifirst) then
            possp=possp+iabs(bytepw)*(ifmin-ifirst)
            ifirst=ifmin
         endif
c
c  Fudge the ending index to prevent use of noise regions.
         ifmax=nint(fmax/delwav)
c         write(*,*)ifmax,ilast,fmax,ilast*delwav
         if(ifmax.lt.ilast) then
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

         zoff=0.0
         zpoff=0.0
         call write_runlog(lunt,col1,specname,iy,doy,gmt,oblat,
     &   oblon,obalt,asza,zpoff,opd,fovi,fovo,amal,ifirst,ilast,
     &   delwav,possp,bytepw,zoff,snr,apf,tins,pins,hins,tout,
     &   pout,hout,lasf,wavtkr,sia,sis,aipl,istat)
         if(mod(ispe,1000).eq.0) write(*,*) ispe
c
      end do ! -------------Main loop over spectra----------------
c
 99   close(luns)
      close(lunt)
      write(*,*) ispe-1, ' spectra found'
      stop
      end
