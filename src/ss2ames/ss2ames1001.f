c  Program SS2AMES converts a space-delimited spreadsheet-format input data
c  file (name.ss) into an Ames-format (1001) output file (name).
c
c  An additional input file (name.hdr) containing further information
c  about the data set (comment lines etc. ) is also read.
c
c  The format of name.ss must be as follows:
c  The first line must contain the number of header lines, and columns.
c  The next header lines will be included as comment lines in the Ames file.
c  The last header line should contain the column labels, one per column.
c
c  Note that the first column of the resulting table (record_#) is inserted 
c  by SS2AMES to be the independent variable in the Ames-format file,
c  guaranteeing that it increase monotonically by a sufficient amount.
c
      integer jcol,kcol,ncol,mcol,ij,j,k,lunss,lunhdr,lunames,
     & nmisscol,nlheadss,ios,ls,lsmax,
     & iobs,nobs,nlheadames,ffi,nscoml,nncoml,lnbc,nvpl
      parameter (
     &   lunss=12,    !  the spreadsheet-format input file
     &   lunhdr=13,   !  the header file
     &   lunames=14,  !  the ames-format output file
     &   mcol=102)    !  maximum number of columns
      integer indxcol(mcol)
      character header*2048,filnam*60,comment*80,string*80,
     & colabel(mcol)*48,misslabel(mcol)*32,version*44,fmt*8,
     & yval(mcol)*12,valmiss(mcol)*12,vscal*12,siobs*6
      parameter (vscal='  1.0000E+00',ffi=1001 )
      filnam=
     & '                                                            '

      version=' SS2AMES   Version 3.0.0    7-Feb-2000   GCT'
      write(6,*)version
      do while (lnbc(filnam) .le. 0)
         write(6,'(a,$)') ' Enter Name of Ames-format Output File : '
         read(5,'(a)') filnam
      end do
c
c  Find out how many columns in input file
      open(lunss,file=filnam(:lnbc(filnam))//'.ss',status='old')
      read(lunss,*) nlheadss,ncol
      nscoml=nlheadss-1
      call skiprec(lunss,nlheadss-3)
      read(lunss,'(8x,a)')string
      call substr(string,valmiss,mcol,ncol)
      read(lunss,'(a)') header
      call substr(header,colabel,mcol,kcol)
      if( kcol .ne. ncol ) stop 'Mismatched number of columns/titles'
c
c  Find out how many rows of data in .ss (input) file
c  and find the longest data-value string (lsmax)
      iobs=0
      ios=0
      lsmax=5
      do while (ios .eq. 0)
         read(lunss,'(a)',iostat=ios) string
         call substr(string,yval,mcol,ncol)
         do jcol=1,ncol
            ls=lnbc(yval(jcol))
            if(ls.gt.lsmax) lsmax=ls
         end do
         iobs=iobs+1
      end do
      nobs=iobs-1
      write(*,*)nobs,lsmax
      nvpl=80/(lsmax+1)  ! Number of Values Per Line
      write(fmt,'(a1,i2,a1,i2.2,a1)') '(',nvpl,'a',lsmax+1,')'
      write(*,*) fmt
c
c  Find out how many comment lines will be needed in output file
      open(lunhdr,file=filnam(:lnbc(filnam))//'.hdr',status='old')
      call skiprec(lunhdr,7)
      read(lunhdr,*) nncoml
      rewind(lunhdr)
c
c  Compute how meany header lines the Ames format Output file will have.
      nlheadames=10+2*(1+(ncol-2)/nvpl)+ncol
     & +(1+nncoml)+(1+nscoml)
c
c  Write Ames-format output file.
      open(lunames,file=filnam,status='unknown')
      write(lunames,*)nlheadames,ffi
      do k=1,7
         read(lunhdr,'(a)') string
         write(lunames,'(a)')string
      end do
      write(lunames,'(a)')' record_#'
      write(lunames,*)ncol
      write(lunames,fmt)(vscal,j=1,ncol)
      write(lunames,fmt)(valmiss(j),j=1,ncol)
      do j=1,ncol
         write(lunames,'(a)') colabel(j)
      end do
c
c  Write the special comment lines from the .ss header.
      rewind (lunss)
      read(lunss,*) nlheadss, ncol
      write(lunames,*)nscoml
      write(lunames,'(a)') ' The following programs created this file:'
      write(lunames,'(a)') version
      do k=3,nscoml
          read(lunss,'(a)') comment
          write(lunames,'(a)') comment
      end do
      call skiprec(lunss,2)  ! skip missing values and header
c
c   Write the normal comment lines from the .hdr file
      read(lunhdr,*) nncoml
      write(lunames,*)nncoml
      do j=1,nncoml
        read(lunhdr,'(a)') string
        write(lunames,'(a)')string
      end do
c
c  Read which gases to cite as missing
      read(lunhdr,'(a)') string
      write(*,'(a)') string
      close(lunhdr)
      call substr(string,misslabel,mcol,nmisscol)
      do k=1,nmisscol
         write(*,*)k,nmisscol,misslabel(k)
         do j=1,ncol
         if( misslabel(k) .eq. colabel(j) ) indxcol(k)=j
         end do
      end do
c
c  Read row of input matrix, insert missing values, and write in Ames format.
      do iobs=1,nobs
        write(siobs,'(i6)') iobs
        read(lunss,'(a)') string
        call substr(string,yval,mcol,ncol)
c
c  Right-shift value
        do j=1,ncol
           ls=lnbc(yval(j))
           yval(j)(lsmax-ls+1:lsmax)=yval(j)(1:ls)
           do k=1,lsmax-ls
             yval(j)(k:k)=' '
           end do
        end do
c
        do j=1,nmisscol
          ij=indxcol(j)
          yval(ij)  = valmiss(ij)     !  Insert missing values into YVAL
          yval(ij+1)= valmiss(ij)     !  Insert missing values into YVAL
        end do
        write(lunames,fmt)siobs,(yval(j),j=1,ncol)
      end do
c
 300  close(lunss)
      close(lunames)
      write(*,*) 'NOBS, NCOL = ',nobs,ncol
      stop
      end

c      include '/ggg/src/comn/substr.f'
c      include '/ggg/src/comn/fbc.f'
c      include '/ggg/src/comn/fnbc.f'
c      include '/ggg/src/comn/lnbc.f'
c      include '/ggg/src/comn/skiprec.f'
