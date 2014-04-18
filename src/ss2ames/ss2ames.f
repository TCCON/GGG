c  Program SS2AMES converts a space-delimited spreadsheet-format input data
c  file (name) into an Ames-format (1001) output file (name.ames).
c
c  An additional input file (name.hdr) containing further information
c  about the data set (comment lines etc. ) is also read.
c
c  The format of name must be as follows:
c  The first line must contain the number of header lines, and columns.
c  The next header lines will be included as comment lines in the Ames file.
c  The last header line should contain the column labels, one per column.
c
c  Note that the first column (record_#) of the resulting table is chosen
c  by SS2AMES to be the independent variable in the Ames-format file,
c  guaranteeing that it increase monotonically by a sufficient amount.
c
      implicit none
      integer jcol,kcol,ncol,mcol,ij,j,k,lunr_ss,lunr_hdr,lunw_ames,
     & nmisscol,nlheadss,ios,ls,lsmax,maxcpl,lf,
     & iobs,nobs,nlheadames,ffi,nscoml,nncoml,lnbc,nvpl
      parameter (
     &   lunr_ss=12,    !  the spreadsheet-format input file
     &   lunr_hdr=13,   !  the header file
     &   lunw_ames=14,  !  the ames-format output file
     &   maxcpl=132,  !  Maximum Characters Per Line
     &   mcol=102)    !  maximum number of columns
      integer indxcol(mcol)
      character header*2048,filnam*256,comment*(maxcpl),string*2048,
     & colabel(mcol)*80,misslabel(mcol)*80,version*64,fmt*8,
     & yval(mcol)*12,valmiss*12,vscal*12,siobs*6
      real*8 ymin,ymax,ymiss,yy
      parameter (vscal='  1.0000E+00', ffi=1001 )

      version=
     &' SS2AMES                  Version 3.21     30-May-2013    GCT  '
      write(6,*)version

      write(filnam,'(80a1)') (' ',j=1,80)
      write(comment,'(200a1)') (' ',j=1,maxcpl)
      do while (lnbc(filnam) .le. 0)
         if (iargc() == 0) then
            write(6,'(a)') ' Enter Name of Input Data File : '
            read(5,'(a)') filnam
         elseif (iargc() == 1) then
            call getarg(1, filnam)
         else
            stop 'Usage: $gggpath/bin/ss2ames amesfile'
         endif
      end do
      lf=lnbc(filnam)
c
c  Find out how many columns in input file
      open(lunr_ss,file=filnam(:lf),status='old')
      read(lunr_ss,*) nlheadss,ncol
      write(*,*)'nlheadss,ncol=',nlheadss,ncol
      nscoml=nlheadss-1
      call skiprec(lunr_ss,nlheadss-4)

c  Missing Values
      read(lunr_ss,'(8x,a)')valmiss   ! missing value
      lsmax=lnbc(valmiss)
      read(valmiss(:lsmax),*)ymiss
c      call substr(string,valmiss,mcol,kcol)
c      if( kcol .ne. ncol ) then
c         write(*,*)'kcol,ncol=',kcol,ncol
c         stop 'wrong number of missing values'
c      endif
c      lsmax=5
c      do jcol=1,ncol
c         ls=lnbc(valmiss(jcol))
c         if(ls.gt.lsmax) lsmax=ls
c      end do

      read(lunr_ss,'(7x,a)') string   ! format
      read(lunr_ss,'(a)') header     ! column headers
      call substr(header,colabel,mcol,kcol)
      if( kcol .ne. ncol ) then
         write(*,*)'kcol,ncol=',jcol,ncol
         stop 'wrong number of column headers'
      endif
c
c  Find out how many rows of data in .ss (input) file
c  and find the longest data-value string (lsmax)
      iobs=0
      ios=0
      do while (ios .eq. 0)
         read(lunr_ss,'(a)',iostat=ios) string
         call substr(string,yval,mcol,ncol)
         do jcol=1,ncol
            ls=lnbc(yval(jcol))
            if(ls.gt.lsmax) lsmax=ls
         end do
         iobs=iobs+1
      end do
      nobs=iobs-1
      write(*,*)'ncol,nobs,lsmax=',ncol,nobs,lsmax
      nvpl=maxcpl/(lsmax+1)  ! Number of Values Per Line
      write(fmt,'(a1,i2,a1,i2.2,a1)') '(',nvpl,'a',lsmax+1,')'
      write(*,*) fmt
c
c  Find out how many comment lines will be needed in output file
      open(lunr_hdr,file=filnam(:lf)//'.hdr',status='old')
      call skiprec(lunr_hdr,8)
      read(lunr_hdr,*) nncoml
      write(*,*)'nncoml=',nncoml
      rewind(lunr_hdr)
c
c  Compute how many header lines the Ames format Output file will have.
      nlheadames=10+2*(1+(ncol-2)/nvpl)+ncol
     & +(1+nncoml)+(1+nscoml)
c
c  Write Ames-format output file.
      open(lunw_ames,file=filnam(:lf)//'.ames',status='unknown')
      read(lunr_hdr,'(a)') string
      write(lunw_ames,'(a)')string(:lnbc(string))
      write(lunw_ames,*)nlheadames,ffi
      do k=1,7
         read(lunr_hdr,'(a)') string
         write(lunw_ames,'(a)')string(:lnbc(string))
      end do
      write(lunw_ames,'(a)')' record_#'
      write(lunw_ames,*)ncol
      write(lunw_ames,fmt)(vscal,j=1,ncol)
      write(lunw_ames,fmt)(valmiss,j=1,ncol)
      do j=1,ncol
         write(lunw_ames,'(a)') colabel(j)(:lnbc(colabel(j)))
      end do
c
c  Write the special comment lines from the .ss header.
      rewind (lunr_ss)
      read(lunr_ss,*) nlheadss, ncol
      write(lunw_ames,*)nscoml
      write(lunw_ames,'(a)')' The following programs created this file:'
      write(lunw_ames,'(a)') version
      do k=3,nscoml
          read(lunr_ss,'(a)') comment
          write(lunw_ames,'(a)') comment(:lnbc(comment))
      end do
      call skiprec(lunr_ss,2)  ! skip missing values and header
c
c   Write the normal comment lines from the .hdr file
      read(lunr_hdr,*) nncoml
      write(*,*)'nncoml=',nncoml
      write(lunw_ames,*)nncoml
      do j=1,nncoml
        read(lunr_hdr,'(a)') string
        write(lunw_ames,'(a)')string(:lnbc(string))
      end do
c
c  Read which gases to cite as missing
      read(lunr_hdr,'(a)') string
      write(*,'(a)') 'missing gases: '//string(:lnbc(string))
      close(lunr_hdr)
      call substr(string,misslabel,mcol,nmisscol)
      do k=1,nmisscol
         write(*,*)k,nmisscol,misslabel(k)
         do j=1,ncol
         if(misslabel(k) .eq. colabel(j)) indxcol(k)=j
         if(misslabel(k)(:lnbc(misslabel(k)))//'_error' .eq. colabel(j))
     &    indxcol(k)=j
         end do
      end do
c
c  Read row of input matrix, insert missing values, and write in Ames format.
      do iobs=1,nobs
        write(siobs,'(i6)') iobs
        read(lunr_ss,'(a)') string
        call substr(string,yval,mcol,ncol)
c
c  Right-shift value
        do j=1,ncol
           ls=lnbc(yval(j))
           yval(j)(lsmax-ls+1:lsmax)=yval(j)(1:ls)
           do k=1,lsmax-ls
             yval(j)(k:k)=' '
           end do
           read(yval(j),*) yy
           if(yy.lt.ymin) ymin=yy
           if(yy.gt.ymax) ymax=yy
        end do
c
        do j=1,nmisscol
          ij=indxcol(j)
          yval(ij)  = valmiss   ! Insert missing value into YVAL
c          yval(ij+1)= valmiss   ! Insert missing value into YVAL_error
        end do
        write(lunw_ames,fmt)siobs,(yval(j),j=1,ncol)
      end do
c
      close(lunr_ss)
      close(lunw_ames)
      write(*,*) 'NOBS, NCOL = ',nobs,ncol
      write(*,*) 'Ymin  Ymax = ',ymin,ymax
      if(ymax.gt.ymiss)
     &  write(*,*)' SS2AMES: Warning: ymax > ymiss',ymax,ymiss
      stop
      end
