c  Program for writing zenith pointing offsets, previously calculated by
c  ZENANG, to the runlog files (both HgCd and InSb).
c
      implicit none
      include "../ggg_int_params.f"

      integer*4 lunr,lunw,lunz,lnbc,lr,istat,ndet,jdet,
     & ifirst,ilast,iyr,iset,possp,bytepw,ispe,ntot,nhl,ncol,i
      parameter (lunr=14,lunw=15,lunz=16)
c
      real*8 year,yearz,aszaz,poutz,totcon,del,
     & totwas,totzen,
     & totwas2,totzen2

      real*8
     & oblat,           ! observation latitude (deg).
     & oblon,           ! observation longitude (deg).
     & obalt,             ! observation altitude (km)
     & zpdtim,           ! Time of ZPD (UT hours)
     & asza,             ! astronomical solar zenith angle (unrefracted)
     & zenoff,           ! zenith pointing offset
     & graw,             ! spacing of raw spectrum (cm-1) 
     & opd,              ! Optical path difference (cm) of interferogram
     & fovi,             ! Internal angular diameter of FOV (radians)
     & fovo,             ! External angular diameter of FOV (radians)
     & amal,             ! angular misalignment of interferometer (radians)
     & zoff,             ! zero-level offset
     & snr,              ! Signal-to-Noise Ratio
     & tins,             ! Temperature INSide the instrument
     & pins,             ! Pressure INSide the instrument
     & hins,             ! Humidity INSide the instrument
     & tout,             ! Temperature OUTside the instrument
     & pout,             ! Pressure OUTside the instrument
     & hout,             ! Humidity OUTside the instrument
     & sia,              ! Solar Intensity (Average)
c    & sis,              ! Solar Intensity (SD)
     & fvsi,             ! Fractional Variation in Solar Intensity
     & wspd,             ! Wind Speed (m/s)
     & wdir,             ! Wind Direction (deg)
     & azim,             ! Solar Azimuth Angle
     & osds,             ! Observer-Sun Doppler Stretch (ppm)
     & aipl,             ! Airmass-Independent Path Length (km)
     & lasf,             ! laser frequency (e.g. 15798.03 cm-1)
     & wavtkr            ! suntracker operating frequency (e.g. 9900 cm-1)

c
      character title*240,
     & col1*1,                    !first column of runlog record
     & apf*2,                     !apodization function (e.g. BX N2, etc)
     & dl*1,
     & gggdir*(mpath),            !ggg directory path (GGGPATH?)
     & specname*(nchar),          !spectrum name
     & rlgfile*120                !name of runlog file

c
      write(6,*) 'ZENCORR    Version 1.3.3     18-Aug-2008     GCT'
      call get_ggg_environment(gggdir, dl)

 1    write(6,'(a)') 'Enter runlog (e.g. flt93avg.brl) '
      read(5,'(a)') rlgfile
      lr=lnbc(rlgfile)
      if(lr.le.0) go to 1
c
c  Open the original runlog and read the first record.
      if(rlgfile(lr-2:lr-2).eq.'b') then
         open(lunr,file=gggdir(:lnbc(gggdir))//'runlogs'//dl//'bal'
     &    //dl//rlgfile, status='old')
         ndet=2
      elseif (rlgfile(lr-2:lr-2).eq.'o') then
         open(lunr,file=gggdir(:lnbc(gggdir))//'runlogs'//dl//'orb'
     &    //dl//rlgfile, status='old')
         ndet=1
      else
         stop 'unknown extention'
      endif
c
c  Open the new runlog
      read(lunr,*) nhl,ncol
      open(lunw,file='new_rl.out',status='unknown')
      write(lunw,*) nhl,ncol
      
      do i=2,nhl
         read(lunr,'(a)') title
         write(lunw,'(a)') title(:lnbc(title))
      end do
c
c  Also open the ZENANG.OUT file containing the pointing offsets
      open(lunz,file='zenang.out',status='old')

      ntot=0
      totwas2=0.0
      totzen2=0.0
      do ispe=1,999999
         read(lunz,*,end=99) yearz,totcon,poutz,aszaz,del
         do jdet=1,ndet
c        call read_runlog(lunr,col1,specname,iyr,iset,zpdtim,
c    &    oblat,oblon,obalt,asza,zenoff,opd,fovi,fovo,amal,ifirst,
c    &    ilast,graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
c    &    tout,pout,hout,lasf,wavtkr,sia,sis,aipl,istat)
         call read_runlog(lunr,col1,specname,iyr,iset,zpdtim,
     &    oblat,oblon,obalt,asza,zenoff,azim,osds,
     &    opd,fovi,fovo,amal,ifirst,ilast,
     &    graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     &    tout,pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)
c         write(*,*)'specname ',specname,totcon,aszaz,asza+zenoff
c         write(*,*)specname,istat,totcon,aszaz,asza+zenoff
          year=iyr+(iset+zpdtim/24)/366.00
         if(abs(abs(year)-abs(yearz)).gt.0.000001) then
            write(*,*)'Year mismatch: ',
     &      specname,year,yearz
         endif
         if(abs(asza+zenoff-aszaz).ge.0.0001) then
            write(*,*)'Zenith angle mismatch: ',
     &      specname,asza+zenoff,aszaz
         endif
         ntot=ntot+1
         totwas=totwas+zenoff
         totwas2=totwas2+zenoff**2
         if(col1.ne.'-') zenoff=zenoff+del
         totzen=totzen+zenoff
         totzen2=totzen2+zenoff**2
c         call write_runlog(lunw,col1,specname,iyr,iset,zpdtim,
c     &   oblat,oblon,obalt,asza,zenoff,opd,fovi,fovo,amal,
c     &   ifirst,ilast,graw,possp,bytepw,zoff,snr,apf,tins,pins,
c     &   hins,tout,pout,hout,lasf,wavtkr,sia,sis,aipl,istat)
         call write_runlog(lunw,col1,specname,iyr,iset,zpdtim,
     &   oblat,oblon,obalt,asza,zenoff,azim,osds,opd,fovi,fovo,amal,
     &   ifirst,ilast,graw,possp,bytepw,zoff,snr,apf,tins,pins,
     &   hins,tout,pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)
         end do  ! jdet=1,ndet
      end do  ! ispe=1,99999
c
 99   close(lunz)
      close(lunr)
      close(lunw)
      write(*,*)'NTOT = ',ntot,ispe-1
      write(*,*)'ZENOFF was = ',totwas/ntot,' +/- ',sqrt(totwas2/ntot)
      write(*,*)'ZENOFF now = ',totzen/ntot,' +/- ',sqrt(totzen2/ntot)
      stop
      end
