c  Program bin2asc.f
c  Converts a series of binary spectra to two-column (freq, signal) ascii format
c
c  Works for any binary spectra having a GGG-format runlog,
c  which is used as the input file.
c
c  Works transparently for I*2 and R*4 binary data types.
c
c  Byte-reversal handled automatically, if the computer
c  that you are working on has a different endian-ness
c  from the computer that wrote the binary spectra.
c
c  ASCII output spectra are created in the local directory.
c
      implicit none
      integer*4
     & lunr,   ! LUN to read input runlogs from
     & luns,   ! LUN to read binary spectras from
     & lunw,   ! LUN to write ascii spectra to
     & mem,    ! maximum buffer size in bytes
     & i,j,
     & iabpw,  ! absolute values of the bytes per word
     & lnbc,   ! function Last Non-Black Character
     & iend,   ! Endianess of host computer
     & npts,   ! Number of spectral values
     & lr      ! 

      parameter (lunr=25,luns=15,lunw=16)
      parameter (mem=4*1024*2048)
      real*4 bufr4(mem/4)
      integer*2  bufi2(mem/2)
      byte bbuf(mem)

      character 
     & runlog*120,inpath*80,chead*1

      integer*4
     & istat,        ! status flag (0=success, 1=EOF)
     & nlhead,       ! 
     & iyr,          ! year
     & iset,         ! day of year
     & ifirst,       ! index of first spectral point in disk file
     & ilast,        ! index of last spectral point in disk file
     & m1, m2,       ! starting and ending spectral point indices in output file
     & possp,        ! Length of attached header in bytes
     & bytepw        ! Bytes per data word, usually 2 (I*2) or 4 (R*4)

      real*8
     & oblat,        ! observation latitude (deg).
     & oblon,        ! observation longitude (deg).
     & obalt,        ! observation altitude (km)
     & asza,         ! astronomical solar/lunar zenith angle (unrefracted)
     & azim,         ! solar/lunar azimuth angle
     & osds,         ! Observer-Sun Doppler Stretch (ppm)
     & opd,          ! Optical path difference (cm) of interferogram
     & graw,         ! spacing of raw spectrum (cm-1) from GETINFO
     & zpdtim,       ! Time of ZPD (UT hours)
     & zenoff,       ! Zenith angle pointing offset (deg)
     & fovi,         ! Internal angular diameter of FOV (radians)
     & fovo,         ! External angular diameter of FOV (radians)
     & amal,         ! angular misalignment of interferometer (radians)
     & zoff,         ! Zero level offset (dimensionless fraction)
     & snr,          ! Signal-to-Noise Ratio (dimensionless)
     & tins,         ! Inside temperature
     & pins,         ! Inside pressure
     & hins,         ! Inside humidity
     & tout,         ! Outside temperature
     & pout,         ! Outside pressure
     & hout,         ! Outside humidity
     & sia,          ! Solar Intensity (Average)
     & fvsi,         ! Fractional Variation of Solar Intensity
     & wspd,wdir,    ! Winf speed & direction
     & aipl,         ! Airmass-Independent Path Length (km)
     & lasf,         ! Laser Frequency (e.g. 15798 cm-1)
     & wavtkr,       ! suntracker frequency (active tracking)
     & nus, nue      ! selected frequency range of interest

      character
     & col1*1,       ! first column of runlog record
     & specname*38,  ! spectrum name
     & version*62,   ! Version number
     & root*80,      ! path to the ggg directories
     & apf*2         ! apodization function (e.g. BX N2, etc)
c
      equivalence (bbuf,bufi2,bufr4)

      write(6,*)
      version=
     &' BIN2ASC                   Version 1.4.3    31-Aug-2009    GCT'
      call getendian(iend)  ! Find endian-ness of host computer

      write(*,*)'Enter path to input file/runlog:'
      read(*,'(a)') runlog
      open(lunr,file=runlog,status='old')
      read(lunr,*) nlhead         ! Skip header line of runlog
      do i=2,nlhead
         read(lunr,*)
      end do

      write(*,*)'Enter Starting & Ending frequencies:'
      write(*,*)'Enter 0 99999 to retain original spectral limits'
      read(*,*) nus,nue

c  Interrogate environmental variable GGGPATH to find location
c  of root partition (e.g. "/home/toon/ggg/" ).
      call getenv('GGGPATH',root) 
      lr=lnbc(root)     ! length of root string (e.g. 14)
c
      istat=0
      do while (istat.eq.0)     ! Main loop over spectra

c  Read input runlog
1        call read_runlog(lunr,col1,specname,iyr,iset,zpdtim,
     &   oblat,oblon,obalt,asza,zenoff,azim,osds,
     &   opd,fovi,fovo,amal,ifirst,ilast,graw,possp,bytepw,zoff,snr,apf,
     &   tins,pins,hins,tout,pout,hout,
     &   sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)
         if(istat.ne.0) exit
 
c  Check that buffer will be large enough
         iabpw=iabs(bytepw)
         m1=int(nus/graw)+1
         if(m1.lt.ifirst) m1=ifirst
         m2=int(nue/graw)
         if(m2.gt.ilast) m2=ilast
         npts=m2-m1+1
         if(npts*iabpw.gt.mem) stop 'Increase parameter NMAX'
         if(npts.lt.1) go to 1

c  Search for binary spectrum "specname"
         call gindfile(root(:lr)//'/config/data_part.lst',specname,
     &   inpath)
         if(lnbc(inpath).eq.0) then
            write(*,*) specname, ' Cant find input spectrum'
            go to 1
         endif
 
c  Open binary spectrum with recl = total length be read
         open(luns,file=inpath,access='direct',status='old',
     &   form='unformatted',recl=possp+iabpw*(m2-ifirst+1))

c  Read spectral header and data values all at once.
         read(luns,rec=1) (chead,j=1,possp+iabpw*(m1-ifirst)),
     &   (bbuf(j),j=1,iabpw*npts)

c  If necessary, byte-reverse data
         if(iend*bytepw.lt.0) call rbyte(bbuf,iabpw,npts)

         close(luns)
  
c  Write ASCI spectrum
         write(6,*)inpath(:lnbc(inpath))
         open(lunw,file='./asc_'//specname,status='unknown')
         write(lunw,*)5,2
         write(lunw,'(a)') version
         write(lunw,'(a)')
     &  '  Spectrum_File_Name                    Year  Day  Hour'//
     &  '   oblat    oblon   obalt    ASZA   POFF    AZIM   OSDS'//
     &  '    OPD   FOVI  FOVO'//
     &  '  AMAL  IFIRST   ILAST    DELTA_NU   POINTER  BPW ZOFF SNR'//
     &  '  APF  tins  pins  hins   tout   pout  hout'//
     &  '  sia  fvsi  wspd  wdir  lasf    wavtkr  aipl'

         call write_runlog(lunw,col1,specname,iyr,iset,zpdtim,oblat,
     &   oblon,obalt,asza,zenoff,azim,osds,opd,fovi,fovo,amal,m1,m2,
     &   graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,tout,
     &   pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)         

         write(lunw,*)' Frequency_(cm-1)  Signal'
         if(iabpw.eq.2) then
            do i=1,npts
              write(lunw,'(f12.6,f9.4)') graw*(i+m1-1),
     &                         float(bufi2(i))/15000.
            end do
         elseif(iabpw.eq.4) then
            do i=1,npts
              write(lunw,'(f12.6,1pe12.4)') graw*(i+m1-1),bufr4(i)
            end do
         else
            stop 'unknown format'
         endif
         close(lunw)

      end do        ! Main loop over spectra
      close(lunr)
      stop
      end
