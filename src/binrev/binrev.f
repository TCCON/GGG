c  Program binrev.f
c
c  Reverses a binary spectrum, e.g., converting a wavelength-sorted file,
c  to a frequency-sorted file, or vice versa.
c  Use to test that GFIT can handle both cases
c
c  Works for any binary spectra having a GGG-format runlog,
c  which is used as the input file. Works for R*4 binary data types.
c
c  Note that the size of the resulting output file can be smaller
c  than that of the input file.  This is because any info following
c  the last spectral point of the input files is discarded.

      implicit none
      integer*4
     & lunr,   ! LUN to read input runlogs from
     & luns,   ! LUN to read binary spectras from
     & mem,    ! maximum buffer size in bytes
     & j,
     & iabpw,  ! absolute values of the bytes per word
     & lnbc,   ! function Last Non-Black Character
     & iend,   ! Endianess of host computer
     & npts,   ! Number of spectral values
     & lr      ! 

      parameter (lunr=14,luns=15)
      parameter (mem=4*1024*2048)
      real*4 bufr4(mem/4)

      character 
     & runlog*120,inpath*80,chead*1

      integer*4
     & istat,        ! status flag (0=success, 1=EOF)
     & iyr,          ! year
     & iset,         ! day of year
     & ifirst,       ! index of first spectral point in disk file
     & ilast,        ! index of last spectral point in disk file
     & possp,        ! Length of attached header in bytes
     & bytepw        ! Bytes per data word, usually 2 (I*2) or 4 (R*4)

      real*8
     & oblat,        ! observation latitude (deg).
     & oblon,        ! observation longitude (deg).
     & obalt,        ! observation altitude (km)
     & asza,         ! astronomical solar zenith angle (unrefracted)
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
     & sis,          ! Solar Intensity (SD)
     & aipl,         ! Airmass-Independent Path Length (km)
     & lasf,         ! Laser Frequency (e.g. 15798 cm-1)
     & wavtkr        ! suntracker frequency (active tracking)

      character
     & col1*1,       ! first column of runlog record
     & runlab*34,    ! spectrum name
     & root*80,      ! path to the ggg directories
     & header*250,   ! 
     & apf*2         ! apodization function (e.g. BX N2, etc)
c
      write(6,*)
     & ' BINREV Program   Version 1.0.0    5-Oct-2008   GCT'
      call getendian(iend)  ! Find endian-ness of host computer

      write(*,*)'Enter path to input file/runlog:'
      read(*,'(a)') runlog
      open(lunr,file=runlog,status='old')
      read (lunr,'(a)') header  

c  Interrogate environmental variable GGGPATH to find location
c  of root partition (e.g. "/home/toon/ggg/" ).
      call getenv('GGGPATH',root) 
      lr=lnbc(root)     ! length of root string (e.g. 14)
c
      istat=0
      do while (istat.eq.0)     ! Main loop over spectra

c  Read input runlog
1        call read_runlog(lunr,col1,runlab,iyr,iset,zpdtim,
     &    oblat,oblon,obalt,asza,zenoff,opd,fovi,fovo,amal,ifirst,
     &    ilast,graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     &    tout,pout,hout,lasf,wavtkr,sia,sis,aipl,istat)
         if(istat.ne.0) exit
 
         iabpw=iabs(bytepw)
         npts=iabs(ilast-ifirst)+1
c
c  Search for binary spectrum "runlab"
         call gindfile(root(:lr)//'/config/data_part.lst',runlab,
     &   inpath)
         if(lnbc(inpath).eq.0) then
            write(*,*) runlab, ' Cant find input spectrum'
            go to 1
         endif
 
c  Read entire spectrum (header and data) all at once.
         open(luns,file=inpath,access='direct',status='old',
     &   form='unformatted',recl=possp+iabpw*npts)
         read(luns,rec=1)(chead,j=1,possp),(bufr4(j),j=1,npts)
         close(luns)
  
c  Write reversed spectrum (header and data) all at once.
         write(6,*)possp+iabpw*npts,inpath(:lnbc(inpath))
         open(luns,file=inpath(:lnbc(inpath))//'.rev',
     &   access='direct',status='unknown',
     &   form='unformatted',recl=possp+iabpw*npts)
         write(luns,rec=1) (chead,j=1,possp),(bufr4(j),j=npts,1,-1)
         close(luns)

      end do        ! Main loop over spectra
      close(lunr)
      stop
      end
