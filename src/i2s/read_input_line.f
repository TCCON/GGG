      subroutine read_input_line(luni,errnum,inpstat,string)
c
c  Input:
c     luni      I*4    Logical Unit Number for the parameter input file
c     errnum    I*4    Program error status: do nothing if not 0 (=ok)
c
c  Input/Output:
c     inpstat   I*4    Value of IOSTAT returned from input file read
c
c  Output:
c     string    C*(*)  String containing the first non-commented line
c
      implicit none

      integer*4
     & luni,       ! Subroutine input argument (see above)
     & errnum,     ! Subroutine input argument (see above)
     & inpstat,    ! Subroutine input/output argument (see above)
     & indexa,     ! General loop index
     & lnbc,ls     ! Integer function Last Non-Blank Character in string

      logical*4
     & replacing   ! True while replacing CTRL-M

      character
     & string*(*)  ! Subroutine output argument (see above)

      if((errnum.eq.0).and.(inpstat.eq.0)) then
         string(1:1)=':'
         do while(((string(1:1).eq.':').or.(lnbc(string).eq.0)).and.
     &          (inpstat.eq.0))
            read(luni,'(a)',iostat=inpstat) string
            ls=lnbc(string)
c
c  Replace CTRL-M from DOS/Win at end of line by a space.
            replacing=( inpstat.eq.0 .and. ls.gt.0 )
            do while(replacing)
               if((ichar(string(ls:ls))).eq.13) then
                  string(ls:ls)=' '
                  ls=lnbc(string)
                  if(ls.eq.0) replacing=.false.
               else
                  replacing=.false.
               endif
            enddo    ! while(replacing)
c
c  If there are in-line comments, remove them
            if( inpstat.eq.0 .and. ls.gt.0 ) then
               indexa=index(string,' :')
               if(indexa.gt.1) then
                  do while(ls.ge.indexa)
                     string(ls:ls)=' '
                     ls=lnbc(string)
                  enddo
               endif
            endif

         enddo      ! while(((string(1:1).eq.':').or.(lnbc(string).eq.0)....
      endif        ! ((errnum.eq.0).and.(inpstat.eq.0))
      return
      end
