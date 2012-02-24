c  Program spectrum_ratioing.f
c
c  Ratios series of binary spectra which are read in from
c  a runlog-format input file (which must include the header lines).
c  Numerator spectra are denoted by a "n", denominator spectra by a "d"
c
c  There can be multiple numerator spectra and/or denominator spectra.
c  Typically, there might be two denominator spectra: acquired before
c  and after the single cell/numerator spectrum. The two denominator
c  spectra will c  be averaged before use. But there is flexibility
c  to handle multiple numerator spectra also.
c
c  The string "ratio" in the input file triggers the ratioing.
c  Prior to that, the numerator and denominator spectra are
c  separately averaged.


c  Outputs ratio spectra in a two-column (freq, signal) ascii format
c  named: asc_rin09xxx.yyy
c  Also creates a series of binary spectra named: rin09xxx.yyy
c  and the associated runlog.
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
      include "../ggg_int_params.f"

      integer*4
     & lunr_rlg,   ! LUN to read input runlogs from
     & lunw_rlg,   ! LUN to write new runlogs to
     & luns,       ! LUN to read binary spectra from
     & lunw_asc,   ! LUN to write ascii spectra to
     & lun_wbs,    ! LUN to write binary spectrum to
     & mem,    ! maximum buffer size in bytes
     & nnum,   ! Number of numerator spectra
     & nden,   ! Number of Denominator spectra
     & i,j,
     & nhl, ncol,
     & iabpw,  ! absolute values of the bytes per word
     & lnbc,   ! function Last Non-Blank Character
     & iend,   ! Endianess of host computer
     & npts,   ! Number of spectral values
     & lr      ! 

      parameter (lunr_rlg=25,lunw_rlg=26,luns=15,lunw_asc=16,lun_wbs=17)
      parameter (mem=4*1024*2048)
      real*4 bufr4(mem/4),buf_num(mem/4),buf_den(mem/4),
     & eps,rat
      integer*2  bufi2(mem/2)
      byte bbuf(mem)

      character 
     & inpath*80,chead*1,
     & col1*1,                    !first column of runlog record
     & apf*2,                     !apodization function (e.g. BX N2, etc)
     & dl*1,
     & gggdir*(mpath),            !ggg directory path (GGGPATH?)
     & specname*(nchar),          !spectrum name
     & version*64,                !current program version
     & rlgfile*120                !name of runlog file

      integer*4
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

      character
     & specname_num*57,  ! spectrum name
     & runlog_header*512, strl*300
c
      equivalence (bbuf,bufi2,bufr4)

      write(6,*)
      version=
     &' spectrum_ratioing         Version 1.0.1    29-Jan-2010    GCT'
      call getendian(iend)  ! Find endian-ness of host computer

      write(*,*)'Enter path to input file/runlog:'
      read(*,'(a)') rlgfile

      eps=0.01
      open(lunr_rlg,file=rlgfile,status='old')
      open(lunw_rlg,file='new_runlog.grl',status='unknown')
      read(lunr_rlg,*) nhl,ncol
      write(lunw_rlg,*) nhl,ncol
      do i=2,nhl-1
         read(lunr_rlg,'(a)') strl
         write(lunw_rlg,'(a)') strl(:lnbc(strl))
      end do
      read(lunr_rlg,'(a)') runlog_header
      write(lunw_rlg,'(a)') runlog_header(:lnbc(runlog_header))

      write(*,*)'Enter Starting & Ending frequencies:'
      write(*,*)'Enter 0 99999 to retain original spectral limits'
      read(*,*) nus,nue

c  Interrogate environmental variable GGGPATH to find location
c  of root partition (e.g. "/home/toon/ggg/" ).
      call get_ggg_environment(gggdir, dl)
      lr=lnbc(gggdir)     ! length of root string (e.g. 14)
c
      nnum=0
      nden=0
      do i=1,mem/4
         buf_num(i)=0.0
         buf_den(i)=0.0
      end do

      istat=0
      do while (istat.eq.0)     ! Main loop over spectra
c  Read input runlog
         read(lunr_rlg,'(a)',end=99) strl
         if (index(strl,'ratio').gt.0) then
c  Write ratio spectrum
            write(6,*) nnum,nden,specname_num
            specname_num(1:1)='r'
            col1=' '
            open(lunw_asc,file='./asc_'//specname_num,status='unknown')
c            write(lunw_asc,*)5,2
c            write(lunw_asc,'(a)') version
c            write(lunw_asc,'(a)') runlog_header
c            write(lunw_asc,'(a)')' Frequency_(cm-1)  Signal'
c           call write_runlog(lunw_asc,col1,specname_num,iyr,iset,zpdtim,
c     &      oblat,oblon,obalt,asza,zenoff,azim,osds,opd,fovi,fovo,amal,
c     &      m1,m2,graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,tout,
c     &      pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)    
           call write_runlog(lunw_rlg,col1,specname_num,iyr,iset,zpdtim,
     &      oblat,oblon,obalt,asza,zenoff,azim,osds,opd,fovi,fovo,amal,
     &      m1,m2,graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,tout,
     &      pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)    

            open(lun_wbs,file=specname_num,access='direct',
     &      status='unknown', form='unformatted',recl=iabpw*npts)

            if(iabpw.eq.2) then
               do i=1,npts
                  rat=(buf_num(i)/nnum+eps)/(buf_den(i)/nden+eps)
                  if(rat.gt.2.0) rat=2.0
                  if(rat.lt.0.0) rat=0.0
                  write(lunw_asc,'(f12.6,f9.5)') graw*(i+m1-1),rat
                  bufi2(i)=nint(15000.0*rat)
               end do
               write(lun_wbs,rec=1) (bufi2(i),i=1,npts)
            elseif(iabpw.eq.4) then
               do i=1,npts
               bufr4(i)=(buf_num(i)/nnum+eps)/(buf_den(i)/nden+eps)
               write(lunw_asc,'(f12.6,1pe12.4)') graw*(i+m1-1),bufr4(i)
               end do
               write(lun_wbs,rec=1) (bufr4(i),i=1,npts)
            else
               stop 'Unrecognized spectral format'
            endif
            close(lunw_asc)
            close(lun_wbs)

            nnum=0
            nden=0
            do i=1,mem/4
               buf_num(i)=0
               buf_den(i)=0
            end do
         else
            backspace(lunr_rlg)
         endif     !  if(index(strl,'ratio').gt.0)
         if(strl(1:1).eq.':') cycle
         call read_runlog(lunr_rlg,col1,specname,iyr,iset,zpdtim,
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
         if(npts.lt.1) cycle

c  Search for binary spectrum "specname"
         call gindfile(gggdir(:lr)//'config'//dl//'data_part.lst',
     &    specname, inpath)
         if(lnbc(inpath).eq.0) then
            write(*,*) specname, ' Cant find input spectrum'
            cycle
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
  
         if(col1.eq.'n') then
            specname_num=specname
            nnum=nnum+1
            if(iabpw.eq.2) then
               do i=1,npts
                  buf_num(i)=buf_num(i)+float(bufi2(i))/15000
               end do
            elseif(iabpw.eq.4) then
               do i=1,npts
                  buf_num(i)=buf_num(i)+bufr4(i)
               end do
            else
               stop 'unknown format'
            endif
            cycle
         elseif(col1.eq.'d') then
            nden=nden+1
            if(iabpw.eq.2) then
               do i=1,npts
                  buf_den(i)=buf_den(i)+float(bufi2(i))/15000
               end do
            elseif(iabpw.eq.4) then
               do i=1,npts
                  buf_den(i)=buf_den(i)+bufr4(i)
               end do
            else
               stop 'unknown format'
            endif
            cycle
         endif

      end do        ! Main loop over spectra
99    close(lunr_rlg)
      close(lunw_rlg)
      stop
      end
