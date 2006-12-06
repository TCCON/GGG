      subroutine fetch(specpath,bytepw,iskip,buf,npts)
c  Reads a contiguous swath of binary data from disk file PATH into array BUF.
c
c INPUTS:
c     PATH    C**  Location of spectrum file.
c     BYTEPW  I*4  Number of bytes per data word (usually +/-2 or +/-4).
c     ISKIP   I*4  Number of data words to be skipped before starting read.
c     NPTS    I*4  Number of data words to be read.
c
c OUTPUTS:
c     BUF(NPTS)     R*4  Array of data words
c
C  Skips the first ISKIP words and then fills BUF with the next NPTS data words.
c  Byte reversing is performed whenever BYTEPW is -ve
c  Data is assumed IEEE binary, except bytepw=5,7,9 which is assumed ASCII.
      implicit none
      INTEGER*4 i4dum,iskip,npts,jpts,bytepw,iend,kk,
     & platform,nlhead,ncol
      character specpath*(*),user*8, rformat*12
      byte bdum(4)
      INTEGER*2 i2dum
      REAL*4 buf(npts),r4dum,freq,tm
c     equivalence (r4dum,i4dum,i2dum,bdum)
      integer*4 iscale    !DG000909
      real*4 xscale
c
c  Determine which kind of computer we are running on; big- or little-endian
      call getendian(iend)      !  iend = +1 on SUN; = -1 on PC
      kk=bytepw*iend            !  kk = +bytepw on Sun; = -bytepw on PC
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c     Platform specification:                     DG000909
c  This code block will hopefully disappear in F95 when the
c  specifications for the unformatted binary read becomes standardized.
      call getenv('LOGNAME',user)
      if(user.ne.'        ')then
         platform=0               !0=Sun, 1=PC-Linux, 2=PC-Win32
         rformat='unformatted'
      else
         platform=2               !0=Sun, 1=PC-Linux, 2=PC-Win32
         user='PC-Win'
         rformat='binary'
      endif
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c
c     Note:                              DG000901
c     Sun and Digital compilers default to form='unformatted' for direct access files
c     Sun compiler functions OK with this default
c     Digital compiler requires form='binary' or values read in are incorrect
c     Format "unformatted" seems to have blocking bytes on PC/Digital Fortran
c     Sun does not accept form='binary'
c
c     bytepw=+/-3 indicates Grams I*4 format spectra
c
      if(iabs(bytepw).eq.+4) then         !  R*4
         open(19,file=specpath,access='direct',status='old',
     &   recl=4,form=rformat)
         xscale=4000000                !Spitzbergen Bruker
         xscale=1                      ! Bruker
         do jpts=1,npts
            read(19,rec=jpts+iskip) r4dum
            if(kk.lt.0) call rbyte(r4dum,4,1)
            buf(jpts)=r4dum/xscale
         end do
c         write(*,*)'fetch: ',kk,buf(1),buf(npts),xscale
c
      elseif(iabs(bytepw).eq.3) then                 !Grams I*4
         open(19,file=specpath,access='direct',status='old',
     &   recl=4,form=rformat)
         read(19,rec=1)bdum
         iscale=bdum(4)
         xscale=float(2**(32-iscale))*50000.
         do jpts=1,npts
            read(19,rec=jpts+iskip) i4dum
            if(kk.lt.0) call rbyte(i4dum,4,1)
            buf(jpts)=float(i4dum)/xscale
         end do
c
      elseif(iabs(bytepw).eq.2) then     !  integer*2
         open(19,file=specpath,access='direct',status='old',
     &   recl=2,form=rformat)
         xscale=20000.
         do jpts=1,npts
            read(19,rec=jpts+iskip) i2dum
            if(kk.lt.0) call rbyte(i2dum,2,1)
            buf(jpts)=float(i2dum)/xscale
         end do
c
      elseif(bytepw.eq.5) then    ! its an ascii file (16i5 format)
         open(19,file=specpath,status='old')
         do jpts=1,iskip/16  ! Skip header & unwanted data
            read(19,*)
         end do
         read(19,'(16f5.0)')(r4dum,jpts=1,mod(iskip,16)), ! skip partial line
     &    (buf(jpts),jpts=1,npts)                    ! then read desired data
         do jpts=1,npts
            buf(jpts)=buf(jpts)/20000.
         end do
c
      elseif(bytepw.eq.7) then    ! its a .spt output file
         open(19,file=specpath,status='old')
         read(19,*) nlhead, ncol
         call skiprec(19,nlhead-1+iskip)  !  Skip header & unwanted data
         do jpts=1,npts  !
            read(19,*)freq,tm,buf(jpts)  ! read wanted data
         end do
c
      elseif(bytepw.eq.9) then    ! its a simple xyplot-format frequency-signal.
         open(19,file=specpath,status='old')
         read(19,*) nlhead, ncol
         call skiprec(19,nlhead-1+iskip)  !  Skip header & unwanted data
         do jpts=1,npts  !
            read(19,*)freq,buf(jpts)  ! read wanted data
         end do
      else
         write(*,*)'BYTEPW=',bytepw
         stop ' FETCH: Unknown file format'
      endif
      close(19)
      RETURN
      END
