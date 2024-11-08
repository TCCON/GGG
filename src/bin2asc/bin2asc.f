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
      include "../gfit/ggg_int_params.f"

      integer*4
     & lunr_rl,   ! LUN to read input runlogs from
     & luns,      ! LUN to read binary spectras from
     & lunw_asc,  ! LUN to write ascii spectra to
     & mem,       ! maximum buffer size in bytes
     & i,j,
     & iabpw,     ! absolute values of the bytes per word
     & lnbc,      ! function Last Non-Blank Character
     & iend,      ! Endianess of host computer
     & npts,      ! Number of spectral values
     & lr         ! 

      parameter (lunr_rl=25,luns=15,lunw_asc=16)
      parameter (mem=4*1024*2048)
      real*4 bufr4(mem/4)
      integer*2  bufi2(mem/2)
      byte bbuf(mem)

      character 
     & fullrlgfile*120,   ! runlog file name with absolute path
     & inpath*190,chead*1,
     & data_fmt_read_rl*256,col_labels_rl*320,
     & cnus*16,cnue*16,
     & col1*1,
     & apf*2,
     & dl*1,
     & gggdir*(mpath),
     & specname*(nchar),
     & version*64

      integer*4 idum,
     & istat,        ! status flag (0=success, 1=EOF)
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

      equivalence (bbuf,bufi2,bufr4)

      idum=mfilepath ! Avoid compiler warning (unused parameter)
      idum=mauxcol ! Avoid compiler warning (unused parameter)
      idum=mcolvav ! Avoid compiler warning (unused parameter)
      idum=mgas    ! Avoid compiler warning (unused parameter)
      idum=mlev    ! Avoid compiler warning (unused parameter)
      idum=mrow_qc ! Avoid compiler warning (unused parameter)
      idum=mspeci  ! Avoid compiler warning (unused parameter)
      idum=mvmode  ! Avoid compiler warning (unused parameter)
      idum=ncell   ! Avoid compiler warning (unused parameter)

      version=
     &' BIN2ASC                Version 1.52        2017-05-22       GCT'
      write(*,*) version
      call getendian(iend)  ! Find endian-ness of host computer

      if (iargc() == 0) then
         write(*,*)'Enter path to input file/runlog:'
         read(*,'(a)') fullrlgfile
         write(*,*)'Enter Starting & Ending frequencies:'
         write(*,*)'Enter 0 99999 to retain original spectral limits'
         read(*,*) nus,nue
      elseif (iargc() == 1) then
         call getarg(1, fullrlgfile)
         nus=0.0d0
         nue=999999.d0
      elseif (iargc() == 3) then
         call getarg(1, fullrlgfile)
         call getarg(2, cnus)
         call getarg(3, cnue)
         read(cnus,*)nus
         read(cnue,*)nue
      else
         stop 'Usage: $gggpath/bin/bin2asc path/runlog [nus, nue]'
      endif

      open(lunr_rl,file=fullrlgfile,status='old')
      call read_runlog_header(lunr_rl,data_fmt_read_rl,col_labels_rl)

c  Interrogate environmental variable GGGPATH to find location
c  of root partition (e.g. "/home/toon/ggg/" ).
      call get_ggg_environment(gggdir, dl)
      lr=lnbc(gggdir)     ! length of root string (e.g. 14)
c
      istat=0
      do while (istat.eq.0)     ! Main loop over spectra

c  Read input runlog
1        call read_runlog_data_record(lunr_rl,data_fmt_read_rl,
     &   col1,specname,iyr,iset,zpdtim,
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
         call gindfile(gggdir(:lr)//'config'//dl//'data_part.lst',
     &     specname, inpath)
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
         write(*,*)inpath(:lnbc(inpath))
         open(lunw_asc,file='./asc_'//specname,status='unknown')
         write(lunw_asc,*)5,2
         write(lunw_asc,'(a)') version(:lnbc(version))
         write(lunw_asc,'(a)') col_labels_rl(:lnbc(col_labels_rl))

         call write_runlog_data_record(lunw_asc,data_fmt_read_rl,
     &   col1,specname,iyr,iset,zpdtim,oblat,
     &   oblon,obalt,asza,zenoff,azim,osds,opd,fovi,fovo,amal,m1,m2,
     &   graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,tout,
     &   pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)         

         write(lunw_asc,*)' Frequency_(cm-1)  Signal'
         if(iabpw.eq.2) then
            do i=1,npts
               write(lunw_asc,'(f12.6,f9.4)') graw*(i+m1-1),
     &         float(bufi2(i))/15000.
            end do
         elseif(iabpw.eq.4) then
            do i=1,npts
               write(lunw_asc,'(f12.6,1pe12.4)') graw*(i+m1-1),bufr4(i)
            end do
         else
            stop 'unknown format'
         endif
         close(lunw_asc)

      end do        ! Main loop over spectra
      close(lunr_rl)
      stop
      end
