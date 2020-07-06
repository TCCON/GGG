      subroutine read_postproc_header
     & (lunr, nlhead, ncol, nrow, nauxcol,
     &  ymiss, data_fmt, addn_lines, naddn)
c
c SUBROUTINE READ_POSTPROC_HEADER
c     Reads in the header of a post-processing file.
c
c   Inputs:
c     lunr          I*4     the file unit to be read from, must be opened
c
c   Outputs:
c     nlhead        I*4     the number of lines in the header
c     ncol          I*4     the total number of columns in the file body
c     nrow          I*4     the number of data rows
c     nauxcol       I*4     the number of columns that contain auxiliary data
c     ymiss         R*8     the value used as a missing value flag
c     data_fmt      C*40    the format of a single row of data in the file body
c     addn_lines    C*(*)   array of additional lines read in from the header
c     naddn         I*4     the number of additional lines read in

      implicit none
      include './postproc_params.f'

      integer*4
     & iline,       ! index of header line
     & lunr,        ! input: file unit
     & nlhead,      ! output: number of header lines
     & ncol,        ! output: number of columns
     & nrow,        ! output: number of non-header rows
     & nauxcol,     ! output: number of auxiliary columns
     & naddn        ! output: number of additional header lines
     
      real*8
     & ymiss        ! output: missing value flag

      character
     & data_fmt*40,                     ! output: format of one data row
     & addn_lines(maddln)*(mcharhead)   ! output: extra lines in the header to propagate forward 


c Start by reading in the number of lines in the header
      read(lunr, countfmt) nlhead,ncol,nrow,nauxcol

c Assume that the last three lines are, in order, the missing value
c flag, the data format, and the column header line. For now we're
c not reading in the column header line in case that needs treated
c specially, but we might be able to change that in the future.
c
c Everything up until the last three lines is considered "additional"
c and is read into the addn_lines array so that it can be written
c out into the next file. This is usually version numbers and 
c correction factors.
      naddn = 0
      do iline=2,nlhead-3
        read(lunr, '(a)') addn_lines(iline-1)
        naddn = naddn + 1
      end do

      read(lunr,'(8x,d12.5)') ymiss
      read(lunr,'(7x,a)') data_fmt

      end subroutine
