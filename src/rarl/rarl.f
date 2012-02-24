c  Program rarl.f (Read ACE RunLog).
c  Reads an ACE-format runlog.
c  Writes a GGG-format runlog.
c  Uses keywords to identify important parameters in ACE runlog (doesn't
c  assume that the order of parameters will always be the same).
c
      implicit none
      include "../ggg_int_params.f"

c  Variables associated with the write-runlog subroutine call.
      integer*4  lsp,lrp,liff,fsib,file_size_in_bytes,
     & lunw,             ! Logical unit number
     & istat,            ! status flag (0=success, 1=EOF)
     & iyr,              ! year
     & iset,             ! day of year
     & ifirst,           ! index of first spectral point in disk file
     & ilast,            ! index of last spectral point in disk file
     & possp,            ! Length of attached header in bytes
     & bytepw            ! Bytes per data word, usually 2 (I*2) or 4 (R*4)

      real*8 ttag,
     & oblat,            ! observation latitude (deg).
     & oblon,            ! observation longitude (deg).
     & obalt,            ! observation altitude (km)
     & tplat,            ! tangent latitude (deg).
     & tplon,            ! tangent longitude (deg).
     & tpalt,            ! tangent altitude (km)
     & asza,l2_asza,     ! astronomical solar zenith angle (unrefracted)
     & opd,              ! Optical path difference (cm) of interferogram
     & graw,             ! spacing of raw spectrum (cm-1) from GETINFO
     & zpdtim,           ! Time of ZPD (UT hours)
     & zenoff,           ! Zenith angle pointing offset (deg)
     & azim,
     & osds,
     & fovi,             ! Internal angular diameter of FOV (radians)
     & fovo,             ! External angular diameter of FOV (radians)
     & amal,             ! angular misalignment of interferometer (radians)
     & zoff,             ! Zero level offset (dimensionless fraction)
     & snr,              ! Signal-To-Noise Ratio (dimensionless)
     & tins,             ! Inside temperature (C)
     & pins,             ! Inside pressure (mbar)
     & hins,             ! Inside humidity (%)
     & tout,             ! Outside temperature (C)
     & pout,             ! Outside pressure (mbar)
     & hout,             ! Outside humidity (%)
     & sia,              ! Solar Intensity (Average)
     & fvsi,             ! Fractional variation in Solar Intensity (SD)
     & wspd,
     & wdir,
     & aipl,             ! Airmass-Independent Path Length (km)
     & lasf,             ! Laser Frequency (e.g. 15798 cm-1)
     & wavtkr            ! Suntracker frequency (active tracking)

      character
     & chiocc*6,
     & col1*1,                    !first column of runlog record
     & apf*2,                     !apodization function (e.g. BX N2, etc)
     & dl*1,
     & gggdir*(mpath),            !ggg directory path (GGGPATH?)
     & specname*(nchar),          !spectrum name
     & version*64                 !current program version

c  Other variables.
      integer*4 lunr,mss,im,id,hh,mm,nhss,ndss,
     & j0,jj,i,k,lnbc,iocc,isun,lrt
      real*8 sec
      parameter (lunr=14,lunw=15,mss=60)
      character header*1600,rlheader*1600,data_string*1600,
     & hss(mss)*32,dss(mss)*32,sp_path*120,
     & ssun(2)*2,rl_path*120,ifformat*6
      logical filexist

      data ssun/'sr','ss'/

      version = ' rarl.f   Version 3.0.2   2011 Aug 22   GCT'
      col1=' '

      call get_ggg_environment(gggdir, dl)
      write(*,*) 'root=',gggdir
      lrt=lnbc(gggdir)
      liff=4
      ifformat='(i4.4)'
      do iocc=5000,38000  ! loop over occultations
      if(iocc.eq.10000) then
        liff=5
        ifformat='(i5.5)'
      endif
      do isun=1,2   !  Sunset/rise
      write(chiocc(1:liff),ifformat)iocc
      rl_path='/export/raid1/spectra/ace/'//ssun(isun)//chiocc(1:liff)//
     & '/runlog.sqltable'
      lrp=lnbc(rl_path)
      inquire(file=rl_path,exist=filexist)
      if(filexist) then
         write(*,*) rl_path(:lrp)//' found'
         open(lunr,file=rl_path,status='old')
         open(lunw, file=gggdir(:lrt)//'runlogs'//dl//
     &   'orb'//dl//ssun(isun)//chiocc(:liff)//'.orl',
     &   status='unknown')
         rlheader=
     &     '     Spectrum_File_Name'//
     &     '                                    Year  Day  Hour'//
     &     '   oblat    oblon   obalt    ASZA    POFF    AZIM   OSDS'//
     &     '    OPD   FOVI  FOVO  AMAL   IFIRST    ILAST    DELTA_NU'//
     &     '   POINTER  BPW ZOFF SNR  APF tins  pins  hins   tout'//
     &     '   pout  hout'//
     &     '  sia    fvsi   wspd  wdir  lasf    wavtkr  aipl'
         call substr(rlheader, hss, mss, nhss)
         if(nhss.gt.mss) stop 'nhss > mss'
         write(lunw,*)3,nhss
         write(lunw,'(a)') version
         write(lunw,'(a)') rlheader(:lnbc(rlheader))

c         write(lunw,'(a)') ' Spectrum_File_Name    Year  Day  Hour'//
c     &     '   oblat    oblon   obalt    ASZA   POFF    OPD   FOVI  FOVO'//
c     &     '  AMAL  IFIRST   ILAST    DELTA_NU   POINTER  BPW ZOFF SNR'//
c     &     '  APF tins  pins  hins   tout   pout  hout   lasf    wavtkr'//
c     &     '  sia   sis   aipl'
         read(lunr,'(a)') header
         call cdsubstr(header,hss,mss,nhss)
         if(nhss.gt.mss) stop 'nhss > mss'
         do i=1,500  ! Loop over spectra in one occultation
            read(lunr,'(a)',end=90) data_string
            call cdsubstr(data_string,dss,mss,ndss)
            l2_asza=0
            do k=1,ndss
               if(hss(k).eq.'timetag')read(dss(k),'(f14.3)',err=88)ttag
               if(hss(k).eq.'obslat') read(dss(k),*,err=88) oblat
               if(hss(k).eq.'obslon') read(dss(k),*,err=88) oblon
               if(hss(k).eq.'obsalt') read(dss(k),*,err=88) obalt
               if(hss(k).eq.'tplat') read(dss(k),*,err=88) tplat
               if(hss(k).eq.'tplon') read(dss(k),*,err=88) tplon
               if(hss(k).eq.'tanht') read(dss(k),*,err=88) tpalt
               if(hss(k).eq.'zenith') read(dss(k),*,err=88) asza
               if(hss(k).eq.'azimuth') read(dss(k),*,err=88) azim
               if(hss(k).eq.'l2_zenith' .and. lnbc(dss(k)) .gt. 0)
     &         read(dss(k),*,err=88) l2_asza
               if(hss(k).eq.'corrected_timestamp') then
                  read(dss(k),'(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f5.2)',
     &            err=88) iyr,im,id,hh,mm,sec
                  call julian(iyr,1,0,j0)
                  call julian(iyr,im,id,jj)
                  iset=jj-j0
                  zpdtim=hh+(dfloat(mm)+sec/60+2.05/60)/60
c add 2.05 seconds to the corrected timestamp to get the ZPD time
               endif
            end do  ! k=1,ndss
c           if(l2_asza .gt. 0) then 
c              zenoff=l2_asza - asza !this takes into account the difference between asza and the ACE-retrieved L2-asza
c           else
               zenoff=0
c           endif
            write(specname,'(a24,f14.3,a4)')
     &      'FTS-TRANSMITTANCE-NA-NA-',ttag,'.adf'
            call gindfile(gggdir(:lrt)//'config'//dl//'data_part.lst',
     &      specname,sp_path)
            lsp=lnbc(sp_path)
            inquire(file=sp_path,exist=filexist)
c            write(*,*) filexist,sp_path
            fsib=file_size_in_bytes(19,sp_path)
            possp=fsib-250001*8
c            write(*,*)'fsib,possp=',fsib,possp
            opd=25.0
            fovi=0.0062
            fovo=0.00125
            amal=0.0058
            ifirst=0
            ilast=250000
            graw=0.02000004883d0
            bytepw=+6
            zoff=0.0
            snr=400.
            apf='BX'
            osds=0.0d0
            tins=28.1
            pins=0.
            hins=99.9
            tout=555.
            pout=0.
            hout=0.
            lasf=6452
            wavtkr=6450
            sia=0.0
            fvsi=0.0
            aipl=0.000
            call write_runlog(lunw,col1,specname,iyr,iset,zpdtim,
     &      oblat,oblon,obalt,asza,zenoff,azim,osds,opd,fovi,fovo,amal,
     &      ifirst,ilast,graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     &      tout,pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)
            if(istat.ne.0) stop 'Error in write_runlog'
88          continue
         end do  ! Loop over spectra
90       close(lunr)
         close(lunw)
      else
c         write(*,*) rl_path(:lrp)//' not found'
      endif  ! if(filexist)
      end do   ! loop over sunset/rise
      end do   ! loop over occultations
      stop
      end
