      subroutine fetch(specpath,bytepw,iskip,buf,npts)

c  Reads a contiguous swath of data from disk file PATH into array BUF.
c
c INPUTS:
c     PATH    C**  Location of spectrum file.
c     BYTEPW  I*4  Number of bytes per data word (usually +/-2 or +/-4).
c     ISKIP   I*4  Number of bytes to be skipped before starting read.
c     NPTS    I*4  Number of data words to be read.
c
c OUTPUTS:
c     BUF(NPTS)     R*4  Array of data words
c
C  Skips the first ISKIP words and then fills BUF with the next NPTS data words.
c  Byte reversing is performed whenever IEND*BYTEPW is -ve
c  Data is assumed IEEE binary, except bytepw=5,7,9 which is assumed ASCII.
      implicit none
      INTEGER*4 i4dum,iskip,npts,j,jpts,bytepw,iend,kk,nlhead,ncol
      character specpath*(*),rformat*12
      byte bdum(4),bb
      INTEGER*2 i2dum
      REAL*4 buf(npts),r4dum,tm
      real*8 freq,fwas,del,dwas
      integer*4 iscale    !DG000909
      real*4 xscale,simag
c
c  Determine which kind of computer we are running on; big- or little-endian
      call getendian(iend)      !  iend = +1 on SUN; = -1 on PC
      kk=bytepw*iend            !  kk = +bytepw on Sun; = -bytepw on PC
C
c     Note:                              DG000901
c     Sun and Digital compilers default to form='unformatted' for direct access files
c     Sun compiler functions OK with this default
c     Digital compiler requires form='binary' or values read in are incorrect
c     Format "unformatted" seems to have blocking bytes on PC/Digital Fortran
c     Sun does not accept form='binary'
      rformat='binary'       ! Windows
      rformat='unformatted'  ! Unix, Linux
c
c     bytepw=+/-3 indicates Grams I*4 format spectra
c
      if(iabs(bytepw).eq.+4) then         !  R*4
c         write(*,*)'iskip,npts=',iskip,npts
         open(19,file=specpath,access='direct',status='old',
     &   recl=iskip+4*npts,form=rformat)
         read(19,rec=1) (bb,j=1,iskip),(buf(j),j=1,npts)
c         write(*,*)'fetch: buf=',(buf(j),j=1,npts)
         if(kk.lt.0) call rbyte(buf,4,npts)
c         do j=1,npts
c           buf(j)=0.001*buf(j)
c         end do
c      if(iabs(bytepw).eq.+4) then         !  R*4
c         open(19,file=specpath,access='direct',status='old',
c     &   recl=4,form=rformat)
c         xscale=4000000                !Spitzbergen Bruker
c         xscale=1                      ! Bruker
c         do jpts=1,npts
c            read(19,rec=jpts+iskip) r4dum
c            if(kk.lt.0) call rbyte(r4dum,4,1)
c            buf(jpts)=r4dum/xscale
c         end do
c
      elseif(iabs(bytepw).eq.3) then                 !Grams I*4
         open(19,file=specpath,access='direct',status='old',
     &   recl=4,form=rformat)
         read(19,rec=1)bdum
         iscale=bdum(4)
         xscale=float(2**(32-iscale))*50000.
         do jpts=1,npts
            read(19,rec=jpts+iskip/4) i4dum
            if(kk.lt.0) call rbyte(i4dum,4,1)
            buf(jpts)=float(i4dum)/xscale
         end do
c
      elseif(iabs(bytepw).eq.2) then     !  integer*2
         open(19,file=specpath,access='direct',status='old',
     &   recl=2,form=rformat)
         xscale=20000.
         do jpts=1,npts
            read(19,rec=jpts+iskip/2) i2dum
            if(kk.lt.0) call rbyte(i2dum,2,1)
            i4dum=i2dum    ! g95
            buf(jpts)=float(i4dum)/xscale
         end do
c
      elseif(bytepw.eq.5) then    ! its an ascii file (16i5 format)
         open(19,file=specpath,status='old')
         do jpts=1,iskip/80  ! Skip header & unwanted data
            read(19,*)
         end do
         read(19,'(16f5.0)')(r4dum,jpts=1,mod(iskip,80)), ! skip partial line
     &    (buf(jpts),jpts=1,npts)                    ! then read desired data
         do jpts=1,npts
            buf(jpts)=buf(jpts)/20000.
         end do
c
c  ACE binary spectra (.adf) containing interleaved R*4 [Real,Imag]
      elseif(iabs(bytepw).eq.+6) then         !  R*4
         open(19,file=specpath,access='direct',status='old',
     &   recl=iskip+8*npts,form=rformat)
         read(19,rec=1) (bb,j=1,iskip),(buf(j),simag,j=1,npts)
         if(kk.lt.0) call rbyte(buf,4,npts)
c
c  The spectral signal is in the third column, e.g. an spt output file
      elseif(bytepw.eq.7) then 
         open(19,file=specpath,status='old')
         read(19,*) nlhead, ncol
c         write(*,*) nlhead, ncol,iskip
         call skiprec(19,nlhead-1+iskip/7)  !  Skip header & unwanted data
         do jpts=1,npts  !
            read(19,*)freq,tm,buf(jpts)  ! read wanted data
c            write(*,*) freq,buf(jpts)
         end do
c
c  The spectral signal is in the second column and there's a header
      elseif(bytepw.eq.9) then 
         fwas=0.0
         dwas=1.0
         open(19,file=specpath,status='old')
         read(19,*) nlhead, ncol
         call skiprec(19,nlhead-1+iskip/bytepw)  !  Skip header & unwanted data
         do jpts=1,npts  !
            read(19,*)freq,buf(jpts)  ! read wanted data
            del=freq-fwas
            if(del*dwas.le.0.0) write(*,*) 'Warning1: FETCH',fwas,freq
            fwas=freq
            dwas=del
         end do
c
c  The spectral signal is in the second column but there's no header
      elseif(bytepw.eq.10) then 
         nlhead=0
         fwas=0.0
         dwas=1.0
         open(19,file=specpath,status='old')
c         read(19,*) nlhead, ncol
         call skiprec(19,nlhead+iskip/bytepw)  !  Skip header & unwanted data
         do jpts=1,npts  !
            read(19,*)freq,buf(jpts)  ! read wanted data
            del=freq-fwas
            if(del*dwas.le.0.0) write(*,*) 'Warning1: FETCH',fwas,freq
            fwas=freq
            dwas=del
c            write(*,*) freq,buf(jpts)
c            if(jpts.eq.1) write(*,*)freq,buf(jpts)
c            if(jpts.eq.npts) write(*,*)freq,buf(npts)
c             buf(jpts)=buf(jpts)*1.e+06
         end do
c  Header, followed by y-values grouped in 4's
      elseif(bytepw.eq.11) then
         open(19,file=specpath,status='old')
         read(19,*) nlhead, ncol
         call skiprec(19,nlhead-1)  !  Skip header
         read(19,*)(r4dum,jpts=1,iskip),(buf(jpts),jpts=1,npts)
      else
         write(*,*)'BYTEPW=',bytepw
         stop ' FETCH: Unknown file format'
      endif
      close(19)
      RETURN
      END
