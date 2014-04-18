      subroutine read_runlog_header(lunr_rl,data_fmt_rl,col_labels_rl)
      integer*4 lunr_rl,mcol,ncol,kcol,nhl,j
      parameter (mcol=40)
      character outarr(mcol)*40,col_labels_rl*(*),data_fmt_rl*(*)
c
c  Inputs:
c    lunr_rl        I*4    ! Logical Unit Number of already-open runlog
c    
c Outputs:
c    data_fmt_rl    C*(*)
c    col_labels_rl  C*(*)
c
c Typical Usage:
c    open(lunw_rlg,file='runlog_name')
c    call write_runlog_header(lunw_rlg,data_fmt_rl,col_labels_rl)
c    do ispec=1,nspec
c       call write)runlog_data_record(lunw_rlg,data_fmt_rl, ....)
c    end do
c    close(lunw_rlg)

c  The format statement below is the default. This allows the code
c  to be able to handle recent runlogs without a format statement.
c  If a format statement is found in the header, it will overwrite
c  the default.  Eventually, when all runlogs have a format statement,
c  the default won't be needed.
      data_fmt_rl='(a1,a57,1x,2i4,f8.4,f8.3,f9.3,2f8.3,1x,f6.4,f8.3,'//
     & 'f7.3,f7.2,3(1x,f5.4),2i9,1x,f14.11,i9,i3,1x,f5.3,i5,1x,'//
     & 'a2,2(f6.1,f8.2,f5.1),f7.1,f7.4,f6.1,f6.0,f10.3,f7.0,f7.3)'

      read(lunr_rl,*)nhl,ncol
      do j=2,nhl
         read(lunr_rl,'(a)') col_labels_rl
         if(col_labels_rl(:7).eq.'format=') then
            data_fmt_rl=col_labels_rl(8:)
         endif
      end do

c The last line of the header is the column labels
      call substr(col_labels_rl, outarr, mcol, kcol)
      if(kcol.ne.ncol) then
         write(*,*)' Error: KCOL, NCOL = ',kcol,ncol
         stop 'read_runlog_header: KCOL .NE. NCOL'
      endif
      return
      end
