      subroutine build_cit_name(outpath,pattern,chnindic,srcindic,
     & year,month,day,ichan,srcode,runno,filename)
c
c  Input:
c    outpath   C*(*)  Directory path to CIT files
c    pattern   C*(*)  Pattern for CIT file-naming convention
c    chnindic  C*(*)  Channel indicators for CIT file-naming convention
c    srcindic  C*(*)  Source indicators for CIT file-naming convention
c    year      I*4    Year from catalog
c    month     I*4    Month from catalog
c    day       I*4    Day from catalog
c    ichan     I*4    Channel number (1=InGaAs=slave, 2=Si=master)
c    srcode    I*4    Numerical code to identify light source selection
c    runno     I*4    Run number
c
c  Output:
c    filename  C*(*)  Full file name for CIT file
c
      implicit none

      integer*4
     & year,       ! Subroutine input argument (see above)
     & month,      ! Subroutine input argument (see above)
     & day,        ! Subroutine input argument (see above)
     & ichan,    ! Subroutine input argument (see above)
     & srcode,     ! Subroutine input argument (see above)
     & runno,      ! Subroutine input argument (see above)
     & lnbc,       ! Integer function Last Non-Blank Character in string
     & pattind,    ! Index into file naming pattern
     & outind,     ! Index into output file name
     & charcnt     ! Count of characters in current sub-pattern

      character
     & outpath*(*), ! Subroutine input argument (see above)       
     & pattern*(*), ! Subroutine input argument (see above)       
     & chnindic*(*),! Subroutine input argument (see above)
     & srcindic*(*),! Subroutine input argument (see above)
     & filename*(*),! Subroutine output argument (see above)
     & currchar,   ! Current character being handled from pattern
     & tempstr*11  ! String used to hold converted numbers

      outind=lnbc(outpath)
      write(filename,'(a)') outpath(1:outind)
      pattind=1

      do while(pattind.le.lnbc(pattern))
         currchar=pattern(pattind:pattind)
         if((currchar.eq.'Y').or.
     &      (currchar.eq.'M').or.
     &      (currchar.eq.'D').or.
     &      (currchar.eq.'R')) then
            charcnt=0
            do while(pattern(pattind:pattind).eq.currchar)
               charcnt=charcnt+1
               pattind=pattind+1
            enddo
            if(charcnt.gt.11) then
               write(*,'(a)')
     &        'Warning: integer format pattern truncated to 11 chars'
               charcnt=11
            endif
            if(currchar.eq.'Y') then
               write(tempstr,'(i11.11)') year
            elseif(currchar.eq.'M') then
               write(tempstr,'(i11.11)') month
            elseif(currchar.eq.'D') then
               write(tempstr,'(i11.11)') day
            elseif(currchar.eq.'R') then
               write(tempstr,'(i11.11)') runno
            endif
            filename(outind+1:outind+charcnt)=tempstr(12-charcnt:11)
            outind=outind+charcnt
         elseif(currchar.eq.'C') then
            outind=outind+1
            filename(outind:outind)=chnindic(ichan:ichan)
            pattind=pattind+1
         elseif(currchar.eq.'S') then
            outind=outind+1
            filename(outind:outind)=srcindic(srcode:srcode)
            pattind=pattind+1
         else
            outind=outind+1
            filename(outind:outind)=currchar
            pattind=pattind+1
         endif
      enddo   ! while(pattind.le.lnbc(pattern))
      return
      end
