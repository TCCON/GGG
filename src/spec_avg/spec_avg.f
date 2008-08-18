c  Program spec_avg.f
c  Averages a series of binary spectra.

c  Works for any binary spectra having a GGG-format runlog,
c  which is used as the input file.
c
c  Works transparently for I*2 and R*4 binary data types.
c
c  Byte-reversal handled automatically, if the computer
c  that you are working on has a different endian-ness
c  from the computer that wrote the binary spectra.
c
c  Headerless binary average output spectra are created
c  in the local directory.
c
      integer*4
     & lun_rpt,   ! LUN to write report
     & lun_rlg,   ! LUN to read input runlogs from
     & luns,   ! LUN to read binary spectras from
     & lunw,   ! LUN to write acerage spectrum to
     & nmax,   ! maximum buffer size in bytes
     & j,
     & iabpw,  ! absolute values of the bytes per word
     & lnbc,   ! function Last Non-Black Character
     & nspe,jspe,i,
     & iend,   ! Endianess of host computer
     & npts,nwas,   ! Number of spectral values
     & lr      ! 

      parameter (lun_rpt=24,lun_rlg=25,luns=15,lunw=16)
      parameter (nmax=4*1024*2048)
      real*4 bufr4(nmax/4),tbuf(nmax/4)
      integer*4 bufi4(nmax/4)
      integer*2  bufi2(nmax/2)
      byte bbuf(nmax)

      character 
     & runlog*120,specpath*80,chead*1

      integer*4
     & istat,        ! status flag (0=success, 1=EOF)
     & iyr,          ! year
     & iset,         ! day of year
     & ifirst,       ! index of first spectral point in disk file
     & ilast,        ! index of last spectral point in disk file
     & m1, m2,       ! starting and ending spectral point indices in output file
     & m1was, m2was, ! starting and ending spectral point indices in output file
     & possp,        ! Length of attached header in bytes
     & bytepw        ! Bytes per data word, usually 2 (I*2) or 4 (R*4)

      real*8 tiyr,tiset,tbytepw,taa,tad,tdd,del,
     & toblat, oblat,        ! observation latitude (deg).
     & toblon, oblon,        ! observation longitude (deg).
     & tobalt, obalt,        ! observation altitude (km)
     & tasza,  asza,         ! astronomical solar zenith angle (unrefracted)
     & topd,   opd,          ! Optical path difference (cm) of interferogram
     & tgraw,  graw,         ! spacing of raw spectrum (cm-1) from GETINFO
     & tzpdtim,zpdtim,       ! Time of ZPD (UT hours)
     & tzenoff,zenoff,       ! Zenith angle pointing offset (deg)
     & tfovi,  fovi,         ! Internal angular diameter of FOV (radians)
     & tfovo,  fovo,         ! External angular diameter of FOV (radians)
     & tamal,  amal,         ! angular misalignment of interferometer (radians)
     & tzoff,  zoff,         ! Zero level offset (dimensionless fraction)
     & tsnr,   snr,          ! Signal-to-Noise Ratio (dimensionless)
     & ttins,  tins,         ! Inside temperature
     & tpins,  pins,         ! Inside pressure
     & thins,  hins,         ! Inside humidity
     & ttout,  tout,         ! Outside temperature
     & tpout,  pout,         ! Outside pressure
     & thout,  hout,         ! Outside humidity
     & tsia,   sia,          ! Solar Intensity (Average)
     & tsis,   sis,          ! Solar Intensity (SD)
     & taipl,  aipl,         ! Airmass-Independent Path Length (km)
     & tlasf,  lasf,         ! Laser Frequency (e.g. 15798 cm-1)
     & twavtkr,wavtkr,       ! suntracker frequency (active tracking)
     & nus, nue      ! selected frequency range of interest

      character
     & col1*1,       ! first column of runlog record
     & runlab*34,    ! spectrum name
     & root*80,      ! path to the ggg directories
     & dplist*80,    ! Data Partition list (~/ggg/config/data_part.lst)
     & apf*2         ! apodization function (e.g. BX N2, etc)
c
      equivalence (bbuf,bufi2,bufi4,bufr4)

      write(6,*)
     & ' SPEC_AVG Program   Version 0.0.1   21-Oct-2007   GCT'

      call getendian(iend)  ! Find endian-ness of host computer

      write(*,*)'Enter path to input_file/runlog:'
      read(*,'(a)') runlog
      open(lun_rlg,file=runlog,status='old')
      read(lun_rlg,*)          ! Skip header line of runlog

      write(*,*)'Enter Starting & Ending frequencies (cm-1):'
      write(*,*)'Enter 0 99999 to retain original spectral limits'
      read(*,*) nus,nue

c  Interrogate environmental variable GGGPATH to find location
c  of root partition (e.g. "/home/toon/ggg/" ).
      call getenv('GGGPATH',root) 
      lr=lnbc(root)     ! length of root string (e.g. 14)
      dplist=root(:lr)//'/config/data_part.lst'
c
c Initialize variables used in the averaging
      jspe=0
      do i=1,nmax/4
         tbuf(i)=0.0
      end do
      tiyr=0.0
      tiset=0.0
      tzpdtim=0.0
      toblat=0.0
      toblon=0.0
      tobalt=0.0
      tasza=0.0
      tzenoff=0.0
      topd=0.0
      tfovi=0.0
      tfovo=0.0
      tamal=0.0
      tgraw=0.0
      tbytepw=0.0
      tzoff=0.0
      tsnr=0.0
      ttins=0.0
      tpins=0.0
      thins=0.0
      ttout=0.0
      tpout=0.0
      thout=0.0
      tlasf=0.0
      twavtkr=0.0
      tsia=0.0
      tsis=0.0
      taipl=0.0

      istat=0
      do while (istat.eq.0)     ! Main loop over spectra

c  Read input runlog
1        call read_runlog(lun_rlg,col1,runlab,iyr,iset,zpdtim,
     &    oblat,oblon,obalt,asza,zenoff,opd,fovi,fovo,amal,ifirst,
     &    ilast,graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     &    tout,pout,hout,lasf,wavtkr,sia,sis,aipl,istat)
         if(istat.ne.0) exit
 
c  Check that buffer will be large enough
         iabpw=iabs(bytepw)
         m1=int(nus/graw)+1
         if(m1.lt.ifirst) m1=ifirst
         m2=int(nue/graw)
         if(m2.gt.ilast) m2=ilast
         npts=m2-m1+1
         if(npts*iabpw.gt.nmax) stop 'Increase parameter NMAX'
         if(npts.lt.1) go to 1

c  Search for binary spectrum "runlab"
         call gindfile(dplist,runlab,specpath)
         if(lnbc(specpath).eq.0) then
            write(*,*) runlab, ' Cant find input spectrum'
            go to 1
         endif
 
c  Check that the spectra contain the same number of points: nus to nue
         if( jspe.gt.0  .and.  npts.ne.nwas ) go to 1
         nwas=npts
         m1was=m1
         m2was=m2

         write(6,*)jspe,specpath(:lnbc(specpath))
c  Open binary spectrum with recl = total length be read
         open(luns,file=specpath,access='direct',status='old',
     &   form='unformatted',recl=possp+iabpw*(m2-ifirst+1))

c  Read spectral header and data values all at once.
         read(luns,rec=1) (chead,j=1,possp+iabpw*(m1-ifirst)),
     &   (bbuf(j),j=1,iabpw*npts)
         close(luns)

c  If necessary, byte-reverse data
         if(iend*bytepw.lt.0) call rbyte(bbuf,iabpw,npts)

         jspe=jspe+1

         do i=1,npts
            tbuf(i)=tbuf(i)+bufr4(i)
         end do
         tiyr= tiyr+ iyr
         tiset= tiset+ iset
         tzpdtim= tzpdtim+ zpdtim
         toblat= toblat+ oblat
         toblon= toblon+ oblon
         tobalt= tobalt+ obalt
         tasza= tasza+ asza
         tzenoff= tzenoff+ zenoff
         topd= topd+ opd
         tfovi= tfovi+ fovi
         tfovo= tfovo+ fovo
         tamal= tamal+ amal
         tgraw= tgraw+ graw
         tbytepw= tbytepw+ bytepw
         tzoff= tzoff+ zoff
         tsnr= tsnr+ snr**2
         ttins= ttins+ tins
         tpins= tpins+ pins
         thins= thins+ hins
         ttout= ttout+ tout
         tpout= tpout+ pout
         thout= thout+ hout
         tlasf= tlasf+ lasf
         twavtkr= twavtkr+ wavtkr
         tsia= tsia+ sia
         tsis= tsis+ sis
         taipl= taipl+ aipl
      end do        ! Main loop over spectra
      close(lun_rlg)
      nspe=jspe

c===================================================
c Divide by nspe to compute average values
      write(*,*)'computing average spectrum'
      do i=1,nwas
         tbuf(i)=tbuf(i)/nspe
      end do
      tiyr=(tiyr+tiset/365.25)/nspe
      iyr= int(tiyr)
      iset= nint(365.25*(tiyr-iyr))
      zpdtim= tzpdtim/nspe
      oblat= toblat/nspe
      oblon= toblon/nspe
      obalt= tobalt/nspe
      asza= tasza/nspe
      zenoff= tzenoff/nspe
      opd = topd/nspe
      fovi= tfovi/nspe
      fovo= tfovo/nspe
      amal= tamal/nspe
      graw= tgraw/nspe
      bytepw= nint(tbytepw/nspe)
      zoff= tzoff/nspe
      snr = sqrt(tsnr)/nspe
      tins= ttins/nspe
      pins= tpins/nspe
      hins= thins/nspe
      tout= ttout/nspe
      pout= tpout/nspe
      hout= thout/nspe
      lasf= tlasf/nspe
      wavtkr= twavtkr/nspe
      sia = tsia/nspe
      sis = tsis/nspe
      aipl= taipl/nspe

c  Write runlog
      write(*,*)'writing average spectrum and runlog'
      open(lunw,file='spec_avg.grl',status='unknown')
      write(lunw,'(a)') ' Spectrum_File_Name    Year  Day  Hour'//
     &  '   oblat    oblon   obalt    ASZA   POFF    OPD   FOVI  FOVO'//
     &  '  AMAL  IFIRST   ILAST    DELTA_NU   POINTER  BPW ZOFF SNR'//
     &  '  APF tins  pins  hins   tout   pout  hout   lasf    wavtkr'//
     &  '  sia   sis   aipl'
      possp=0
      runlab='spec_avg.out'
      call write_runlog(lunw,col1,runlab,iyr,iset,zpdtim,oblat,
     &oblon,obalt,asza,zenoff,opd,fovi,fovo,amal,m1was,m2was,
     &graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,tout,
     &pout,hout,lasf,wavtkr,sia,sis,aipl,istat)         
      close(lunw)

c  If necessary, byte-reverse data
      if(iend*bytepw.lt.0) call rbyte(tbuf,4,nwas)
c
c  Write headerless binary spectrum
      write(*,*)nwas
      open(lunw,file='spec_avg.out',access='direct',status='unknown',
     &form='unformatted',recl=4*nwas)
      write(lunw,rec=1) (tbuf(j),j=1,nwas)
      close(lunw)

c==================================================================
c  Check how closely the individual spectra compare with the average.
      write(*,*)'Comparing individual spectra with average...'
      open(lun_rpt,file='spec_avg.rpt',status='unknown')
      open(lun_rlg,file=runlog,status='old')
      read(lun_rlg,*)          ! Skip header line of runlog
      do jspe=1,nspe     !  Loop over spectra

c  Read input runlog
2        call read_runlog(lun_rlg,col1,runlab,iyr,iset,zpdtim,
     &    oblat,oblon,obalt,asza,zenoff,opd,fovi,fovo,amal,ifirst,
     &    ilast,graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     &    tout,pout,hout,lasf,wavtkr,sia,sis,aipl,istat)
         if(istat.ne.0) exit ! end of runlog
 
c  Check that buffer will be large enough
         iabpw=iabs(bytepw)
         m1=int(nus/graw)+1
         if(m1.lt.ifirst) m1=ifirst
         m2=int(nue/graw)
         if(m2.gt.ilast) m2=ilast
         npts=m2-m1+1
         if(npts.lt.1) go to 2

c  Search for binary spectrum "runlab"
         call gindfile(dplist,runlab,specpath)
         if(lnbc(specpath).eq.0) then
            write(*,*) runlab, ' Cant find input spectrum'
            go to 2
         endif
 
c  Check that the spectra contain the same number of points: nus to nue
         if( jspe.gt.0  .and.  npts.ne.nwas ) go to 2
         nwas=npts
         m1was=m1
         m2was=m2

c         write(6,*)specpath(:lnbc(specpath))
c  Open binary spectrum with recl = total length be read
         open(luns,file=specpath,access='direct',status='old',
     &   form='unformatted',recl=possp+iabpw*(m2-ifirst+1))

c  Read spectral header and data values all at once.
         read(luns,rec=1) (chead,j=1,possp+iabpw*(m1-ifirst)),
     &   (bbuf(j),j=1,iabpw*npts)
         close(luns)

         taa=0.0
         tad=0.0
         tdd=0.0
         do i=1,npts
            del=bufr4(i)-tbuf(i)
            taa=taa+bufr4(i)**2
            tad=tad+bufr4(i)*del
            tdd=tdd+del**2
         end do
         write(lun_rpt,'(a,2f9.4)')runlab,1.+tad/taa,
     &   sqrt(tdd*taa-tad*tad)/taa

      end do ! jspe=1,nspe     ! Loop over spectra
      close(lun_rlg)
      close(lun_rpt)

      stop
      end
