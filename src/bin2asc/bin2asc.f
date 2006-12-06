c  Program bin2asc.f
c
c  Converts binary spectra to ascii format.
c  Works for any binary spectra having a GGG-format runlog.
c  Uses the appropriate GGG-format runlog as an input file.
c  Works transparently for I*2 and R*4 data types.
c  Byte-reversal handled automatically, if needed.
c
c   To run, type, e.g.
c       /home/toon/ggg/bin/bin2asc < pa20050309.grl
c  ASCII output spectras are created in the local directory
c
      integer*4 lunr,luns,lunw,nmax,j,iabpw,lnbc,i,iend,spsv,ispe,lr
c     &,kk,irec
      parameter (lunr=5,luns=15,lunw=16)
      parameter (nmax=4*1024*2048)
c
      integer*4 buf4(nmax/4)
c
      integer*2  buf2(nmax/2)
c
      byte bbuf(nmax)

      character 
     & inpath*32,chead*1

      integer*4
     & istat,            ! status flag (0=success, 1=EOF)
     & iyr,              ! year
     & iset,             ! day of year
     & ifirst,           ! index of first spectral point in disk file
     & ilast,            ! index of last spectral point in disk file
     & possp,            ! Length of attached header in bytes
     & bytepw            ! Bytes per data word, usually 2 (I*2) or 4 (R*4)

      real*8
     & oblat,            ! observation latitude (deg).
     & oblon,            ! observation longitude (deg).
     & obalt,             ! observation altitude (km)
     & asza,             ! astronomical solar zenith angle (unrefracted)
     & opd,              ! Optical path difference (cm) of interferogram
     & graw,             ! spacing of raw spectrum (cm-1) from GETINFO
     & zpdtim,           ! Time of ZPD (UT hours)
     & zenoff,           ! Zenith angle pointing offset (deg)
     & fovi,             ! Internal angular diameter of FOV (radians)
     & fovo,             ! External angular diameter of FOV (radians)
     & amal,             ! angular misalignment of interferometer (radians)
     & zoff,             ! Zero level offset (dimensionless fraction)
     & snr,              ! Signal-to-Noise Ratio (dimensionless)
     & tins,             ! Inside temperature
     & pins,             ! Inside pressure
     & hins,             ! Inside humidity
     & tout,             ! Outside temperature
     & pout,             ! Outside pressure
     & hout,             ! Outside humidity
     & sia,              ! Solar Intensity (Average)
     & sis,              ! Solar Intensity (SD)
     & aipl,             ! Airmass-Independent Path Length (km)
     & lasf,             ! Laser Frequency (e.g. 15798 cm-1)
     & wavtkr            ! suntracker frequency (active tracking)

      character
     & col1*1,           ! first column of runlog record
     & runlab*21,        ! spectrum name
     & root*64,          ! path to the ggg directories
     & apf*2             ! apodization function (e.g. BX N2, etc)
c
      equivalence (bbuf,buf2,buf4)

      write(6,*)
     & ' BIN2ASC Program   Version 1.2.1   11-Oct-2006   GCT'

      call getendian(iend)  ! Find endian-ness of host computer
      read(lunr,*)          ! Skip header line of runlog
c
      do ispe=1,999999  ! Main loop over spectra

         call read_runlog(lunr,col1,runlab,iyr,iset,zpdtim,
     &    oblat,oblon,obalt,asza,zenoff,opd,fovi,fovo,amal,ifirst,
     &    ilast,graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     &    tout,pout,hout,lasf,wavtkr,sia,sis,aipl,istat)
         if(istat.ne.0) exit
 
         call getenv('GGGPATH',root)
         lr=lnbc(root)
         call gindfile(root(:lr)//'/config/data_part.lst',runlab,
     &   inpath)
         if(lnbc(inpath).eq.0) then
            write(*,*) runlab
            stop ' Cant find input spectrum'
         endif
 
         iabpw=iabs(bytepw)
         spsv=ilast-ifirst+1
         if(iabpw*spsv.gt.nmax) stop 'iabpw*spsv>nmax'
c  Read spectral header and data values all at once.
         open(luns,file=inpath,access='direct',status='old',
     &   form='unformatted',recl=possp+iabpw*spsv)
         read(luns,rec=1) (chead,j=1,possp),(bbuf(j),j=1,iabpw*spsv)
         if(iend*bytepw.lt.0) call rbyte(bbuf,iabpw,spsv)
         close(luns)
  
         write(6,*)ispe,'  '//inpath
         open(lunw,file='./asc_'//runlab,status='unknown')
         write(lunw,*)3,2
         write(lunw,*)runlab,zpdtim,oblat,oblon,obalt,pout,tout
         write(lunw,*)' Frequency_(cm-1)  Signal'
         if(iabpw.eq.2) then
            do i=1,spsv
               write(lunw,'(f9.4,i7)') graw*(i+ifirst-1),buf2(i)
            end do
         elseif(iabpw.eq.4) then
            do i=1,spsv
               write(lunw,'(f11.5,1pe12.4)') graw*(i+ifirst-1),buf4(i)
            end do
         else
            stop 'unknown format'
         endif
         close(lunw)
      end do   ! loop over spectra
      close(lunr)
      stop
      end
