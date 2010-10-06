         subroutine rravgcom(mode,nus,nue,lun_rlg,lun_wav,lun_rpt,
     &   ntype,jtype,nspe,krec,istat)

c  Reads the runlog-format input file from the current position down
c  to the next "average" marker.
c  Mode 0: Does nothing functional (testing only). Simply reads runlog. 
c  Mode 1: Averages the spectrum and header values, writes averages. 
c  Mode 2: Compares individual spectra with the average, writing
c          intensity scale factors and their rms deviations.

      implicit none

c----------------------------------------------------------
c  The following are the MkIV spectral header parameters
c
      character blank*1640,comments*70,detctr*4,detsn*4,filter*4,
     & ippd*12,ippnam*6,preamp*4,rund*2,runinfo*18,runloc*4,sigchn*4,
     & strt*12,com*70
c
      integer*2 fftpow,ialias,idecim,irun,iset,iyr,phase,sfct
c
      integer*4 pinl,pspl,spsv,sspp,totp
c
      real*4 alat,alon,altd,ccor,fovr,gains(8),hins,hout,ippver,
     & offsts(8),pins,pinv,pott,pout,psym,rzero,sampl,soaz,soze,
     & spwn,stnr,tins,tout,wdir,wspd,zmst,zpdtim,zpdv,zsym
c
      real*8 lasf,pspv,sins,sinv,ssps,sspv,zpdl
c----------------------------------------------------------

      integer*4
     & mode,        ! What to do with the spectra: average or compare
     & lun_rlg,     ! LUN to read input runlogs from
     & lun_rbs,     ! LUN to read binary spectra from
     & lun_wbs,     ! LUN to write binary spectra to
     & lun_wav,     ! LUN to write average runlog
     & lun_rpt,     ! lun_rpt+jtype = LUN of file to write rms deviations
     & nmax,        ! maximum buffer size in bytes
     & iabpw,       ! absolute values of the bytes per word
     & lnbc,        ! function Last Non-Black Character
     & ntype,jtype, ! Periodicity of spectra to be averaged
     & krec,        ! number of records read (including invalid ones)
     & nspe,        ! Number of spectra in current average
     & iend,        ! Endianess of host computer
     & kpts,npts,   ! Number of spectral values
     & la,ls,lr,    ! string lengths
     & header_length_mkiv, ! MkIV header length
     & hedlen,
     & i,j      ! 

      parameter (lun_rbs=19,lun_wbs=20,nmax=4*1024*2048,
     & header_length_mkiv=2048)
      real*4 bufr4(nmax/4),tbuf(nmax/4),yi
      integer*4 bufi4(nmax/4)
      integer*2  bufi2(nmax/2)
      byte bbuf(nmax)

      character 
     & cdum*1,
     & specpath*120,strl*512

      integer*4
     & istat,        ! status flag (0=success, 1=EOF)
     & irr,          ! status flag (0=success, 1=EOF)
     & iyrrl,        ! year (from runlog)
     & idoy,         ! day of year
     & ifirst,       ! index of first spectral point in disk file
     & ilast,        ! index of last spectral point in disk file
     & k1, k2,       ! starting and ending spectral point indices in output file
     & n1, n2, ! starting and ending spectral point indices in output file
     & possp,        ! Length of attached header in bytes
     & bytepw        ! Bytes per data word, usually 2 (I*2) or 4 (R*4)

      real*8 tiyrrl,tidoy,tbytepw,taa,tad,tdd,del,pi,d2r,
     & ttotp, taltd, tpout, ttout, thout, tpott, tsoze, tsfct,
     & tspwn, tpspv, tfovr, tzpdtim, talat, talon, tstnr, trzero,
     & toblat, oblat,        ! observation latitude (deg).
     & toblon, oblon,        ! observation longitude (deg)
     & tobalt, obalt,        ! observation altitude (km)
     & tasza,  asza,         ! astronomical solar zenith angle (unrefracted)
     & tazim,  azim,         ! solar azimuth angle (deg)
     & tosds,  osds,         ! observer-sun doppler stretch (ppm)
     & topd,   opd,          ! Optical path difference (cm) of interferogram
     & tgraw,  graw,         ! spacing of raw spectrum (cm-1) from GETINFO
     & tzpdtimrl,zpdtimrl,   ! Time of ZPD (UT hours)
     & tzenoff,zenoff,       ! Zenith angle pointing offset (deg)
     & tfovi,  fovi,         ! Internal angular diameter of FOV (radians)
     & tfovo,  fovo,         ! External angular diameter of FOV (radians)
     & tamal,  amal,         ! angular misalignment of interferometer (radians)
     & tzoff,  zoff,         ! Zero level offset (dimensionless fraction)
     & tsnr,   snr,          ! Signal-to-Noise Ratio (dimensionless)
     & tr8tins,  r8tins,         ! Inside temperature
     & tr8pins,  r8pins,         ! Inside pressure
     & tr8hins,  r8hins,         ! Inside humidity
     & tr8tout,  r8tout,         ! Outside temperature
     & tr8pout,  r8pout,         ! Outside pressure
     & tr8hout,  r8hout,         ! Outside humidity
     & tsia,   sia,          ! Solar Intensity (Average)
     & tfvsi,  fvsi,         ! Solar Intensity (SD)
     & r8wspd,         ! Wind Speed (m/s)
     & r8wdir,         ! Wind Direction (deg.)
     & twe,    twn,          ! Wind vectors (East & North)
     & taipl,  aipl,         ! Airmass-Independent Path Length (km)
     & tlaserf,  laserf,         ! Laser Frequency (e.g. 15798 cm-1)
     & twavtkr,wavtkr,       ! suntracker frequency (active tracking)
     & nus, nue      ! selected frequency range of interest

      character
     & col1*1,       ! first column of runlog record
     & specname*38,  ! spectrum name
     & avgspecname*38,  ! average spectrum name
     & last_valid*38,  ! last valid spectrum name
     & root*80,      ! path to the ggg directories
     & dplist*80,    ! Data Partition list (~/ggg/config/data_part.lst)
     & apf*2         ! apodization function (e.g. BX N2, etc)

      equivalence (bbuf,bufi2,bufi4,bufr4)

      save tbuf  ! save average spectrum between mode=1 and mode=2 calls

      call getendian(iend)  ! Find endian-ness of host computer

c  Interrogate environmental variable GGGPATH to find location
c  of root partition (e.g. "/home/toon/ggg/" ).
      call getenv('GGGPATH',root)
      lr=lnbc(root)     ! length of root string (e.g. 14)
      dplist=root(:lr)//'/config/data_part.lst'

      pi=4*datan(1.0d0)
      d2r=pi/180

      if (mode.eq.1) then
         do j=1,nmax/4
            tbuf(j)=0.0
         end do
         tiyrrl=0.0
         tidoy=0.0
         tzpdtimrl=0.0
         toblat=0.0
         toblon=0.0
         tobalt=0.0
         tasza=0.0    
         tosds=0.0    
         tzenoff=0.0
         tazim=0.0
         topd=0.0
         tfovi=0.0    
         tfovo=0.0
         tamal=0.0    
         tgraw=0.0
         tbytepw=0.0
         tzoff=0.0
         tsnr=0.0
         tr8tins=0.0
         tr8pins=0.0
         tr8hins=0.0
         tr8tout=0.0 
         tr8pout=0.0 
         tr8hout=0.0
         tlaserf=0.0 
         twavtkr=0.0
         tsia=0.0
         tfvsi=0.0
         twe=0.0
         twn=0.0
         taipl=0.0
c
c  Zero MkIV header parameters
        ttotp=0
        taltd=0.0
        tpout=0.0
        ttout=0.0
        thout=0.0
        tpott=0.0
        tsoze=0.0
        tsfct=0.0
        tspwn=0.0
        tpspv=0.0
        tfovr=0.0
        tzpdtim=0.0
        talat=0.0
        talon=0.0
        tstnr=0.0
        trzero=0.0

      endif ! mode=1
c
      n1=0
      n2=0
      krec=0     ! Number of runlog reads
      nspe=0      ! number of valid spectra read
      do         ! Loop over single spectra
c         write(*,*) 'mode, krec, nspe=', mode,krec,nspe
            read(lun_rlg,'(a)',end=99) strl
c         irec=irec+1
c         write(*,*)'strl=',strl(:22)
         krec=krec+1
         if (index(strl,'average').gt.0) exit 
         if(mod(krec-1,ntype).eq.jtype-1) then  ! the right flavor detector
            if(strl(1:1).eq.':') cycle
            backspace(lun_rlg)
            call read_runlog(lun_rlg,col1,specname,
     &      iyrrl,idoy,zpdtimrl,oblat,oblon,obalt,asza,
     &      zenoff,azim,osds,opd,fovi,fovo,amal,
     &      ifirst,ilast,graw,possp,bytepw,zoff,snr,apf,
     &      r8tins,r8pins,r8hins,r8tout,r8pout,r8hout,sia,fvsi,
     &      r8wspd,r8wdir,laserf,wavtkr,aipl,irr)
            if(irr.ne.0) stop 'Error reading runlog (wrong format)'
            if(col1.eq.':') cycle
            if(nspe.eq.0) avgspecname=specname

Check that buffer will be large enough
            iabpw=iabs(bytepw)
            k1=int(nus/graw)+1
            if(k1.lt.ifirst) k1=ifirst
            k2=int(nue/graw)
            if(k2.gt.ilast) k2=ilast
            kpts=k2-k1+1
            if(kpts*iabpw.gt.nmax) stop 'Increase parameter NMAX'
            if(kpts.lt.1) cycle

c  Search for binary spectrum "specname"
            call gindfile(dplist,specname,specpath)
            if(lnbc(specpath).eq.0) then
               write(*,*)' Cant find input spectrum: '//specname
               cycle
            endif

c  Skip spectra containing different starting or ending indices from the first.
            if( nspe.eq.0 ) then
               n1=k1
               n2=k2
               npts=kpts
            else
               if (k1.ne.n1 .or. k2.ne.n2) cycle
            endif

            last_valid=specname 

c  Open binary spectrum with recl = total length be read
            open(lun_rbs,file=specpath,access='direct',status='old',
     &      form='unformatted',recl=possp+iabpw*(k2-ifirst+1))

c  Read spectral header and data values all at once.
            if(iabpw.eq.2) then
              hedlen=header_length_mkiv
              read(lun_rbs,rec=1) runinfo,strt,runloc,iset,irun,rund,
     &  ialias,pins,
     &  tins,hins,lasf,altd,alat,alon,soaz,soze,tout,pout,hout,wspd,
     &  wdir,pinv,pinl,psym,zpdl,zpdv,zsym,zmst,sampl,totp,detctr,
     &  filter,preamp,sigchn,fovr,sinv,sins,ippd,fftpow,idecim,sspp,
     &  spsv,spwn,stnr,sspv,ssps,phase,pspl,pspv,sfct,ippver,comments,
     &  iyr,ccor,pott,zpdtim,detsn,ippnam,rzero,gains,offsts,blank, 
     &        (cdum,j=header_length_mkiv+1,possp+iabpw*(k1-ifirst)),
     &        (bbuf(j),j=1,iabpw*npts)
            else
              hedlen=0
              read(lun_rbs,rec=1) (cdum,j=1,possp+iabpw*(k1-ifirst)),
     &        (bbuf(j),j=1,iabpw*npts)
            endif
            close(lun_rbs)

c  If necessary, byte-reverse data
            if(iend*bytepw.lt.0) call rbyte(bbuf,iabpw,npts)

            ls=lnbc(specname)
            if( mode.eq.2 ) then   ! find rms deviation from average
               if(nspe.eq.0)
     &         write(*,*)'Computing RMS deviations from mean'
               taa=0.0
               tad=0.0
               tdd=0.0
               do i=1,npts
                  if(iabpw.eq.2) then
                    yi=bufi2(i)
                  else
                    yi=bufr4(i)
                  endif
                  del=yi-tbuf(i)
                  taa=taa+yi**2
                  tad=tad+yi*del
                  tdd=tdd+del**2
               end do
               write(lun_rpt+jtype,'(a,2f9.4)')specname(:ls),
     &         1.+tad/taa,sqrt(tdd*taa-tad*tad)/taa
               write(*,'(a,2f9.4)') specpath(:lnbc(specpath)),
     &         1.+tad/taa,sqrt(tdd*taa-tad*tad)/taa
            elseif (mode.eq.1) then  ! Compute average spectrum and runlog
c               write(*,*)'rravgcom: ',mode,jtype,krec,' '//specname(:ls)
               write(*,'(a)') specpath(:lnbc(specpath))

               if(iabpw.eq.2) then
                  do i=1,npts
                     tbuf(i)=tbuf(i)+bufi2(i)
                  end do
               elseif(iabpw.eq.4) then
                  do i=1,npts
                     tbuf(i)=tbuf(i)+bufr4(i)
                  end do
               else
                  write(*,*) 'rravgcom: unsupported data format'
               endif
               tiyrrl= tiyrrl+ iyrrl
               tidoy= tidoy+ idoy
               tzpdtimrl= tzpdtimrl+ zpdtimrl
               toblat= toblat+ oblat
               toblon= toblon+ oblon
               tobalt= tobalt+ obalt
               tasza= tasza + asza
               tazim= tazim + azim
               tosds= tosds + osds
               tzenoff= tzenoff+ zenoff
               topd= topd+ opd
               tfovi= tfovi+ fovi
               tfovo= tfovo+ fovo
               tamal= tamal+ amal
               tgraw= tgraw+ graw
               tbytepw= tbytepw+ bytepw
               tzoff= tzoff+ zoff
               tsnr= tsnr+ snr**2
               tr8tins= tr8tins+ r8tins
               tr8pins= tr8pins+ r8pins
               tr8hins= tr8hins+ r8hins
               tr8tout= tr8tout+ r8tout
               tr8pout= tr8pout+ r8pout
               tr8hout= tr8hout+ r8hout
               tlaserf= tlaserf+ laserf
               twavtkr= twavtkr+ wavtkr
               tsia= tsia+ sia
               tfvsi= tfvsi+ fvsi
               twe=twe+r8wspd*dsin(d2r*r8wdir)
               twn=twn+r8wspd*dcos(d2r*r8wdir)
               taipl= taipl+ aipl
c MkIV header parameters
               ttotp=ttotp+totp
               taltd=taltd+altd
               tpout=tpout+pout
               ttout=ttout+tout
               thout=thout+hout
               tpott=tpott+pott
               tsoze=tsoze+soze
               tsfct=tsfct+sfct
               tspwn=tspwn+spwn
               tpspv=tpspv+pspv
               tfovr=tfovr+fovr
               tzpdtim=tzpdtim+zpdtim
               talat=talat+alat
               talon=talon+alon
               tstnr=tstnr+stnr**2
               trzero=trzero+rzero

            endif ! mode.eq.1
            nspe=nspe+1
         endif  !  if(mod(krec,ntype).eq.jtype) then
      end do   ! Loop over single spectra


      if(mode.eq.1 .and. nspe.gt.0) then
         la=lnbc(avgspecname)
         avgspecname=avgspecname(:la)//'_'//last_valid(la-2:la)
         write(*,*)'Computing average of these ',nspe,' spectra'
         write(lun_rpt+jtype,*) '    average: '//avgspecname
c Divide by nspe to compute average values
         do i=1,npts
            tbuf(i)=tbuf(i)/nspe
         end do
         tiyrrl=(tiyrrl+tidoy/365.25)/nspe
         iyrrl= int(tiyrrl)
         idoy= nint(365.25*(tiyrrl-iyrrl))
         zpdtimrl= tzpdtimrl/nspe
         oblat= toblat/nspe
         oblon= toblon/nspe
         obalt= tobalt/nspe
         asza= tasza/nspe
         azim= tazim/nspe
         osds= tosds/nspe
         zenoff= tzenoff/nspe
         opd = topd/nspe
         fovi= tfovi/nspe
         fovo= tfovo/nspe
         amal= tamal/nspe
         graw= tgraw/nspe
         bytepw= nint(tbytepw/nspe)
         zoff= tzoff/nspe
         snr = sqrt(tsnr)/nspe
         r8tins= tr8tins/nspe
         r8pins= tr8pins/nspe
         r8hins= tr8hins/nspe
         r8tout= tr8tout/nspe
         r8pout= tr8pout/nspe
         r8hout= tr8hout/nspe
         laserf= tlaserf/nspe
         wavtkr= twavtkr/nspe
         sia = tsia/nspe
         fvsi = tfvsi/nspe
         r8wspd = sqrt(twe**2+twn**2)/nspe
         r8wdir = datan2(twe,twn)/d2r
         if(r8wdir.lt.0) r8wdir=r8wdir+360
         aipl= taipl/nspe
         iabpw=iabs(bytepw)

c  MkIV header parameters
        totp=ttotp/nspe
        altd=taltd/nspe
        pout=tpout/nspe
        tout=ttout/nspe
        hout=thout/nspe
        pott=tpott/nspe
        soze=tsoze/nspe
        sfct=tsfct/nspe
        spwn=tspwn/nspe
        pspv=tpspv/nspe
        fovr=tfovr/nspe
        zpdtim=tzpdtim/nspe
        alat=talat/nspe
        alon=talon/nspe
        stnr=sqrt(tstnr)
        rzero=trzero/nspe

         call write_runlog(lun_wav,col1,avgspecname,iyrrl,idoy,zpdtimrl,
     &   oblat,oblon,obalt,asza,zenoff,azim,osds,opd,fovi,fovo,amal,
     &   n1,n2,graw,hedlen,iabpw*iend,zoff,snr,apf,
     &   r8tins,r8pins,r8hins,r8tout,r8pout,r8hout,
     &   sia,fvsi,r8wspd,r8wdir,laserf,wavtkr,aipl,istat)

cc  If necessary, byte-reverse data
c         if(iend*bytepw.lt.0) then call rbyte(tbuf,4,npts)
c
c  Write headerless binary average spectrum
         write(*,*)'writing average spectrum: '//
     &   avgspecname(:lnbc(avgspecname))
         open(lun_wbs,file=avgspecname,access='direct',
     &   status='unknown',form='unformatted',recl=hedlen+iabpw*npts)
         if(iabpw.eq.4) then
            write(lun_wbs,rec=1) (tbuf(j),j=1,npts)
         else
            do i=1,npts
              bufi2(i)=nint(tbuf(i))
            end do
            sspp=n1
            spsv=npts
            ippver=6.0    ! version number
            rund='A '
            ippnam='sp_avg6'
            com='avg of ## : '//specname
            write(com(8:9),'(i2.2)') nspe
            comments=com
            write(lun_wbs,rec=1) runinfo,strt,runloc,iset,irun,rund,
     &  ialias,pins,
     &  tins,hins,lasf,altd,alat,alon,soaz,soze,tout,pout,hout,wspd,
     &  wdir,pinv,pinl,psym,zpdl,zpdv,zsym,zmst,sampl,totp,detctr,
     &  filter,preamp,sigchn,fovr,sinv,sins,ippd,fftpow,idecim,sspp,
     &  spsv,spwn,stnr,sspv,ssps,phase,pspl,pspv,sfct,ippver,comments,
     &  iyr,ccor,pott,zpdtim,detsn,ippnam,rzero,gains,offsts,blank,
     &  (bufi2(j),j=1,npts)
         endif
         close(lun_wbs)

      endif   ! nspe.gt.0
 
      return
99    istat=1
      write(*,*)' rravgcom: Hit EOF of input file'
      return
      end
