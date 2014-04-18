c  Compare Opus Binary Spectra
c  Program compares the spectrum part of pairs of OPUS files and 
c  calculates the rms differences. The header is skipped.
      implicit none
      integer*4 lunr_in,lunw,lunr_bs,possp(2),msp,nsp(2),
     & j,dtype,lp(2),iver,lnblnk,
     & lr,kspec,nspec
      parameter (lunr_in=12,lunr_bs=14,lunw=15,msp=960000,dtype=1031)
      real*4 r4buf(msp,2)
      real*8 tdif,tsum,tt
      character datapart(2)*128,infile*80,
     & specname*80,specpath*256,cdum*1

c  Get name of input file 
      lr=0
      do while (lr .le. 0)
         if (iargc() == 0) then
            write(6,'(a)') 'Input File (e.g. roh.in): '
            read(*,'(a)')infile
         elseif (iargc() == 1) then
            call getarg(1, infile)
         else
            stop 'Use: $gggpath/bin/roh infile '
         endif
         lr=lnblnk(infile)
      end do   !  while (lr .le. 0)

c  Open input file and read the two spectrum paths.
      open(lunr_in,file=infile,status='old')
      read(lunr_in,'(a)') datapart(1)
      read(lunr_in,'(a)') datapart(2)
      lp(1)=lnblnk(datapart(1))
      lp(2)=lnblnk(datapart(2))

      open(lunw,file=infile(:lr)//'.out',status='unknown')
      write(lunw,*)'4 7'
      write(lunw,'(a)') datapart(1)(:lp(1))
      write(lunw,'(a)') datapart(2)(:lp(2))
      write(lunw,*)'    Spectrum           POSSP1   POSSP2   NSP1   NSP2
     &     NOISE      DIFF'

c  Main Loop: reading the names of the spectra
      do kspec=1,999999  ! Loop over spectra
         read(lunr_in,'(a)',end=99) specname

         do iver=1,2   ! Loop over versions/paths

c   Form full specrum path
            specpath=datapart(iver)(:lp(iver))//specname

c   Read OPUS header index to find location & length of spectrum
c   inside OPUS file
            call read_opus_hi(lunr_bs,specpath,dtype,
     &      possp(iver),nsp(iver))

c  Read full binary spectrum
            open(lunr_bs,file=specpath,access='direct',status='old',
     &       form='unformatted',recl=possp(iver)+4*nsp(iver))
            read(lunr_bs,rec=1) (cdum,j=1,possp(iver)),
     &      (r4buf(j,iver),j=1,nsp(iver))
            close(lunr_bs)

         end do  ! Loop over versions/paths

c  Compute mean spectral differences (within a spectrum and
c  between different spectra), normalized by mean spectral values
         tdif=0.0d0
         tsum=0.0d0
         tt=0.0d0
         do j=2,min(nsp(1),nsp(2))-1
            tt=tt+abs(r4buf(j,1)-0.5*(r4buf(j+1,1)+r4buf(j-1,1)))
            tdif=tdif+abs(r4buf(j,2)-r4buf(j,1))
            tsum=tsum+abs(r4buf(j,1))
         end do

c  Write results to output file.
         write(lunw,'(a22,4i8,2(1pe12.4))') specname(:22),
     &   possp(1),possp(2),nsp(1),nsp(2),
     &   0.7071*tt/tsum,tdif/tsum

      end do  ! kspec=1,999999
99    nspec=kspec-1
c     write(*,*)'Number of spectra=',nspec
      stop
      end

      subroutine read_opus_hi(lunr_bs,path,dtype,ipoint,ilen)
c  Reads the OPUS file header index to find the position
c  (byte offset) and length of the binary spectrum.
c
c  Inputs:
c     LUNS  I*4  Logical Unit Number for Binary Spectrum
c     PATH  C**  Full path to spectrum
c    DTYPE  I*4  Data type for spectrum (1031)
c
c  Outputs:
c    IPOINT  I*4  Byte offset
c    ILEN    I*4  Number of R*4 spectral points
c
      integer*4 lunr_bs,mdb,ndb,ip,itype,dtype,ilen,ipoint,
     & reclen,magic
      parameter (reclen=12)
      real*8 prog
      character*120 path

c  Read Header Block.
      open(lunr_bs,file=path,form='unformatted',status='old',
     & access='direct',recl=reclen)

      read(lunr_bs,rec=1)magic,prog      ! magic number, version number
c      write(*,*)magic,prog

      read(lunr_bs,rec=2)ip,mdb,ndb      ! pointer, max. size, cur.  size
c      write(*,*)'ip, mdb, ndb =', ip,mdb,ndb

      ip=ip/reclen  ! convert from bytes to records
      do i=1,ndb
         read(lunr_bs,rec=i+ip) itype,ilen,ipoint
         if(mod(itype,32768).eq.dtype) go to 88   ! Spectrum block found
      end do

      ipoint=0  ! Byte offset
      ilen=0    ! Number of data points
88    close(lunr_bs)
      return
      end
