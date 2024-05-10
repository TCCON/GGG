      subroutine write_postproc_header
     & (lunw, ncol, nrow, nauxcol, ymiss,
     &  data_fmt, addn_lines, naddn, nextra)
c
c SUBROUTINE WRITE_POSTPROC_HEADER
c   Writes the standard header for a post processing output file. Note
c   that the number of header lines is NOT required as an input because
c   it can be inferred
c
c   Inputs:
c       lunw            I*4     the file unit to write to, must be opened
c       ncol            I*4     the number of data columns in the file body
c       nrow            I*4     the number of data rows
c       nauxcol         I*4     the number of columns that contain
c                                 auxiliary data
c       ymiss           R*8     the value used as a missing value flag
c       data_fmt        C*40    the format of a single row of data in
c                                 the file body
c       addn_lines      C*(*)   array of extra lines in the header.
c                                 These will go immediately after the
c                                 first line with the number of header
c                                 lines, etc.
c       naddn           I*4     how many additional lines there are to
c                                 write in addn_lines
c       nextra          I*4     how many additional header lines will be
c                                 written outside of this function. Will
c                                 be added to the nhead number.
      implicit none
      include './postproc_params.f'

      integer*4
     & iline,       ! index of header line
     & lunw,        ! input: file unit
     & nlhead,      ! input: number of header lines
     & ncol,        ! input: number of columns
     & nrow,        ! input: number of non-header rows
     & nauxcol,     ! input: number of auxiliary columns
     & naddn,       ! input: number of additional header lines in addn_lines
     & nextra       ! input: number of extra header lines to be manually written
     
      real*8
     & ymiss        ! input: missing value flag

      character
     & data_fmt*40,                     ! input: format of one data row
     & addn_lines(maddln)*(mcharhead)   ! input: extra lines in the header to propagate forward 

c     The total number of lines in the header is the additional
c     lines plus the first one, the missing value, the data format,
c     and the column headers (4 minimum)
      nlhead = naddn + 4 + nextra
      write(lunw, countfmt) nlhead, ncol, nrow, nauxcol
      do iline=1,naddn
        write(lunw, '(a)') trim(addn_lines(iline))
      end do
      write(lunw, '(a8,1pe12.4)') 'missing:', ymiss
      write(lunw, '(a7,a40)') 'format:', data_fmt

      end subroutine
