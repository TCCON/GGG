c  Program reverse_spectrum (formerly binrev.f)
c
c  Converts a wavelength-sorted spectrum into a frequency-sorted
c  spectrum, or vice versa. Used to test that the GGG code produces
c  consistent results for both file types.
c
c  Works for any binary spectra having a GGG-format runlog,
c  which is used as the input file. Works for R*4 binary data types.
c
c  Note that the size of the resulting output file can be smaller
c  than that of the input file.  This is because any info following
c  the last spectral point of the input file is lost.
c
c  This program was written to check that GFIT correctly handles
c  wavelength-sorted spectra, in additional to frequency-sorted
c  spectra. Of course, once this is verified (by reversing a
c  series of spectra and confirming that the frequency- and 
c  wavelength-sorted versions give identic results), then this
c  program immediately becomes obsolete. But it is good to
c  occasionally check that nothing bad has happened to the code
c  to undo its ability to correctly handle reversed spectra.
c
c   Input Files:
c       runlog (user prompted for path)
c       spectra (found using gindfile)
c
c  Output files:
c       new runlog (in local directory)
c       new spectra (in same directory as originals)
c
c  New spectra are named: oldspectrum.rev
c

      implicit none
      include "../gfit/ggg_int_params.f"

      integer*4
     & lunr_rl,   ! LUN to read input runlogs from
     & lunw_rl,   ! LUN to write output runlogs to
     & lunr_spec, ! LUN to read binary spectras from
     & lunw_spec, ! LUN to write reversed binary spectra to
     & mem,       ! maximum buffer size in bytes
c     & nlhead,ncol,
     & j,idum,
     & iabpw,     ! absolute values of the bytes per word
     & lnbc,      ! function Last Non-Black Character
     & iend,      ! Endianess of host computer
     & npts,      ! Number of spectral values
     & lr         ! 

      parameter (lunr_rl=14,lunw_rl=15,lunr_spec=16,lunw_spec=17)
      parameter (mem=8*1024*2048)
      real*4 bufr4(mem/4)

      character 
     & inpath*80,chead*1,
     & col1*1,                    !first column of runlog record
     & apf*2,                     !apodization function (e.g. BX N2, etc)
     & dl*1,
     & version*56,
     & data_fmt_read_rl*256,
     & data_fmt_write_rl*256,
     & col_labels_rl*320,
     & gggdir*(mpath),            !ggg directory path (GGGPATH?)
     & specname*(nchar),          !spectrum name
     & rlgfile*120                !name of runlog file


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
     & asza,         ! astronomical solar zenith angle (deg)
     & azim,         ! solar azimuth angle (deg)
     & osds,         ! Observer-Sun Doppler Stretch
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
     & fvsi,         ! Fractional Variation of Solar Intensity ()
     & wdir,         ! Wind Direction (deg)
     & wspd,         ! Wind Speed (m/s)
     & aipl,         ! Airmass-Independent Path Length (km)
     & lasf,         ! Laser Frequency (e.g. 15798 cm-1)
     & wavtkr        ! suntracker frequency (active tracking)

      idum=mfilepath ! Avoid compiler warning (unused parameter)
      idum=mauxcol  ! Avoid compiler warning (unused parameter)
      idum=mcolvav  ! Avoid compiler warning (unused parameter)
      idum=mgas     ! Avoid compiler warning (unused parameter)
      idum=mlev     ! Avoid compiler warning (unused parameter)
      idum=mrow_qc  ! Avoid compiler warning (unused parameter)
      idum=mspeci   ! Avoid compiler warning (unused parameter)
      idum=mvmode   ! Avoid compiler warning (unused parameter)
      idum=ncell    ! Avoid compiler warning (unused parameter)
      idum=nchar    ! Avoid compiler warning (unused parameter)

      version=
     & ' REVERSE_SPECTRUM     Version 1.11    2016-04-02    GCT '
      write(6,*)version
      call getendian(iend)  ! Find endian-ness of host computer

      write(*,*)'Enter path to runlog:'
      read(*,'(a)') rlgfile
      open(lunr_rl,file=rlgfile,status='old')
      open(lunw_rl,file='binrev.out',status='unknown')
      call read_runlog_header(lunr_rl,data_fmt_read_rl,col_labels_rl)
c      read(lunr_rl,*) nlhead,ncol
c      do j=2,nlhead
c         read (lunr_rl,*)
c      end do
      call write_runlog_header(lunw_rl,version,data_fmt_write_rl)

c  Interrogate environmental variable GGGPATH to find location
c  of root partition (e.g. "/home/toon/ggg/" ).
      call get_ggg_environment(gggdir, dl) 
      lr=lnbc(gggdir)     ! length of root string (e.g. 14)
c
      istat=0
      do while (istat.eq.0)     ! Main loop over spectra

c  Read input runlog
1       call read_runlog_data_record(lunr_rl,data_fmt_read_rl,
     &  col1,specname,iyr,iset,zpdtim,
     &  oblat,oblon,obalt,asza,zenoff,azim,osds,
     &  opd,fovi,fovo,amal,ifirst,ilast,graw,possp,bytepw,zoff,snr,apf,
     &  tins,pins,hins,tout,pout,hout,
     &  sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)

         if(istat.ne.0) exit
         iabpw=iabs(bytepw)
         npts=iabs(ilast-ifirst)+1
c
c  Search for binary spectrum "specname"
         call gindfile(gggdir(:lr)//'config'//dl//'data_part.lst',
     &     specname,inpath)
         if(lnbc(inpath).eq.0) then
            write(*,*) specname, ' Cant find input spectrum '//specname
            go to 1
         endif
 
c  Write new runlog record. 
        idum=ifirst
        ifirst=-ilast
        ilast=-idum
        graw=-graw
        specname=specname(:lnbc(specname))//'.rev'
        call write_runlog_data_record(lunw_rl,data_fmt_write_rl,
     &  col1,specname,iyr,iset,zpdtim,
     &  oblat,oblon,obalt,asza,zenoff,azim,osds,
     &  opd,fovi,fovo,amal,ifirst,ilast,graw,possp,bytepw,zoff,
     &  snr,apf,tins,pins,hins,tout,pout,hout,
     &  sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)

c  Read entire spectrum (header and data) all at once.
         open(lunr_spec,file=inpath,access='direct',status='old',
     &   form='unformatted',recl=possp+iabpw*npts)
         read(lunr_spec,rec=1)(chead,j=1,possp),(bufr4(j),j=1,npts)
         close(lunr_spec)
  
c  Write reversed spectrum (header and data) all at once.
         write(6,*)possp+iabpw*npts,inpath(:lnbc(inpath))
         open(lunw_spec,file=inpath(:lnbc(inpath))//'.rev',
     &   access='direct',status='unknown',
     &   form='unformatted',recl=possp+iabpw*npts)
         write(lunw_spec,rec=1) (chead,j=1,possp),(bufr4(j),j=npts,1,-1)
         close(lunw_spec)

      end do        ! Main loop over spectra
      close(lunr_rl)
      close(lunw_rl)
      stop
      end
