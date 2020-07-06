      subroutine get_opus_i4(luns,ipoint,param,iend,path,verbose,
     & slicnt,errnum,val)
c
c  This subroutine will search in the file opened at 'luns' the parameter
c  block at byte offset 'ipoint' for the integer*4 named in 'param'.  The
c  returned value will be in 'val' and any detected errors will be flagged
c  in 'errnum'.
c
c  This subroutine checks that parameters of subsequent slices match those
c  of the first slice: this is why it needs 'slicnt' as its input argument.
c  To avoid this check, simply pass the value 1 instead of 'slicnt'.
c
c  Argument 'iend' is used to trigger a byte swap on big-endian machines.
c  Arguments 'path' and 'verbose' are used only for stdout messages.
c
c  Input:
c    luns      I*4    Logical Unit Number for slice input file
c    ipoint    I*4    Pointer to parameter block (in byte units)
c    param     C*(3)  Name of parameter to search for
c    iend      I*4    Computer endianness flag
c    path      C*(*)  Path to IFS125 interferogram slice file
c    verbose   I*4    Level of verbosity for displayed messages
c    slicnt    I*4    Slice counter, used for parameter consistency checks
c
c  Output:
c    errnum    I*4    Error code (0=ok, <0=fatal, >0=recoverable)
c
c  Input/Output:
c    val       I*4    Value of parameter read from OPUS header
c
      implicit none

      integer*4
     & luns,       ! Subroutine input argument (see above)
     & ipoint,     ! Subroutine input argument (see above)
     & iend,       ! Subroutine input argument (see above)
     & verbose,    ! Subroutine input argument (see above)
     & slicnt,     ! Subroutine input argument (see above)
     & errnum,     ! Subroutine output argument (see above)
     & val,        ! Subroutine input/output argument (see above)
     & lnbc,       ! Integer function Last Non-Blank Character in string
     & i4val       ! Value of parameter, interpreted as INT32

      integer*2
     & ers         ! Expected size of the parameter value, in 16-bit units
      parameter (ers=2)

      integer*2
     & i2val(ers), ! Value of parameter, interpreted as array of INT16
     & rs,         ! Returned size of the parameter value, in 16-bit units
     & typ         ! Type of the parameter, read from file

      character
     & param*3,    ! Subroutine input argument (see above)
     & path*(*)    ! Subroutine input argument (see above)

      equivalence (i2val,i4val)

      call get_opus_prm(luns,ipoint,param,iend,ers,i2val,rs,typ)
      if(rs.ne.ers) then
         errnum=4
         write(*,'(4a)')'Warning: ',
     &   path(1:lnbc(path)),' has no value for ',param
      elseif(typ.ne.0) then
         errnum=4
         write(*,'(4a)')'Error: ',
     &   path(1:lnbc(path)),' has bad type for ',param
c      elseif(slicnt.eq.1) then
      elseif(slicnt.ge.1) then
         val=i4val
         if(verbose.ge.4) then
            write(*,*) param,'=',val
         endif
      elseif(val.ne.i4val) then
         errnum=4
         write(*,'(4a)')'Error: get_opus_i4:',
     &   path(1:lnbc(path)),' has inconsistent value for ',param
      endif

      return
      end
c ===================================================================
      subroutine get_opus_r8(luns,ipoint,param,iend,path,verbose,
     & slicnt,errnum,val)
c
c  This subroutine will search in the file opened at 'luns' the parameter
c  block at byte offset 'ipoint' for the real*8 named in 'param'.  The
c  returned value will be in 'val' and any detected errors will be flagged
c  in 'errnum'.
c
c  This subroutine checks that parameters of subsequent slices match those
c  of the first slice: this is why it needs 'slicnt' as its input argument.
c  To avoid this check, simply pass the value 1 instead of 'slicnt'.
c
c  Argument 'iend' is used to trigger a byte swap on big-endian machines.
c  Arguments 'path' and 'verbose' are used only for stdout messages.
c
c  Input:
c    luns      I*4    Logical Unit Number for slice input file
c    ipoint    I*4    Pointer to parameter block (in byte units)
c    param     C*(3)  Name of parameter to search for
c    iend      I*4    Computer endianness flag
c    path      C*(*)  Path to IFS125 interferogram slice file
c    verbose   I*4    Level of verbosity for displayed messages
c    slicnt    I*4    Slice counter, used for parameter consistency checks
c
c  Output:
c    errnum    I*4    Error code (0=ok, <0=fatal, >0=recoverable)
c
c  Input/Output:
c    val       R*8    Value of parameter read from OPUS header
c
      implicit none

      integer*4
     & luns,       ! Subroutine input argument (see above)
     & ipoint,     ! Subroutine input argument (see above)
     & iend,       ! Subroutine input argument (see above)
     & verbose,    ! Subroutine input argument (see above)
     & slicnt,     ! Subroutine input argument (see above)
     & errnum,     ! Subroutine output argument (see above)
     & lnbc        ! Integer function Last Non-Blank Character in string

      integer*2
     & ers         ! Expected size of the parameter value, in 16-bit units
      parameter (ers=4)

      integer*2
     & i2val(ers), ! Value of parameter, interpreted as array of INT16
     & rs,         ! Returned size of the parameter value, in 16-bit units
     & typ         ! Type of the parameter, read from file

      real*8
     & val,        ! Subroutine input/output argument (see above)
     & r8val       ! Value of parameter, interpreted as REAL64

      character
     & param*3,    ! Subroutine input argument (see above)
     & path*(*)    ! Subroutine input argument (see above)

      equivalence (i2val,r8val)

      call get_opus_prm(luns,ipoint,param,iend,ers,i2val,rs,typ)
      if(rs.ne.ers) then
         errnum=4
         write(*,'(4a)')'Warning: ',
     &   path(1:lnbc(path)),' has no value for ',param
      elseif(typ.ne.1) then
         errnum=4
         write(*,'(4a)')'Error: ',
     &   path(1:lnbc(path)),' has bad type for ',param
      elseif(slicnt.ge.1) then
         val=r8val
         if(verbose.ge.4) then
            write(*,*) param,'=',val
         endif
      elseif(val.ne.r8val) then
         errnum=4
         write(*,*) ' slicnt, val, r8val = ',slicnt, val,r8val
         write(*,'(4a)')'Error: get_opus_r8: ',
     &   path(1:lnbc(path)),' has inconsistent value for ',param
      endif

      return
      end
c ===================================================================
      subroutine get_opus_string(luns,ipoint,param,exptyp,minlen,iend,
     & path,verbose,slicnt,errnum,length,cval)
c
c  This subroutine will search in the file opened at 'luns' the parameter
c  block at byte offset 'ipoint' for the character string named in 'param'.
c  The returned null-terminated string will be in 'cval' and its length,
c  excluding the terminator, in 'length'.  Any detected errors will be
c  flagged in 'errnum'.
c
c  Argument 'minlen' is used for the length check: if this minimum
c  value is not met, an error message is produced and 'errnum' is set.
c
c  This subroutine checks that parameters of subsequent slices match those
c  of the first slice: this is why it needs 'slicnt' as its input argument.
c  To avoid this check, simply pass the value 1 instead of 'slicnt'.
c
c  Argument 'iend' is used to trigger a byte swap on big-endian machines.
c  Arguments 'path' and 'verbose' are used only for stdout messages.
c
c  Input:
c    luns      I*4    Logical Unit Number for slice input file
c    ipoint    I*4    Pointer to parameter block (in byte units)
c    param     C*(3)  Name of parameter to search for
c    exptyp    I*4    Expected parameter type (string/enum/senum)
c    minlen    I*4    Minimum expected length for string answer (excl. NULL)
c    iend      I*4    Computer endianness flag
c    path      C*(*)  Path to IFS125 interferogram slice file
c    verbose   I*4    Level of verbosity for displayed messages
c    slicnt    I*4    Slice counter, used for parameter consistency checks
c
c  Output:
c    errnum    I*4    Error code (0=ok, <0=fatal, >0=recoverable)
c    length    I*4    Length of string returned in cval, excluding NULL term.
c
c  Input/Output:
c    cval      C*(*)  NULL-terminated value of param read from OPUS header
c
      implicit none

      integer*4
     & luns,       ! Subroutine input argument (see above)
     & ipoint,     ! Subroutine input argument (see above)
     & exptyp,     ! Subroutine input argument (see above)
     & minlen,     ! Subroutine input argument (see above)
     & iend,       ! Subroutine input argument (see above)
     & verbose,    ! Subroutine input argument (see above)
     & slicnt,     ! Subroutine input argument (see above)
     & errnum,     ! Subroutine output argument (see above)
     & length,     ! Subroutine output argument (see above)
     & lnbc,       ! Integer function Last Non-Blank Character in string
     & indexa,     ! General loop index
     & indexb      ! General loop index

      integer*2
     & mrs         ! Maximum size of the parameter value, in 16-bit units
c      parameter (mrs=24)
      parameter (mrs=48)  ! GCT 2014-06-15

      integer*2
     & i2val(mrs), ! Value of parameter, interpreted as array of INT16
     & rs,         ! Returned size of the parameter value, in 16-bit units
     & typ         ! Type of the parameter, read from file

      integer*1
     & bval(2*mrs) ! Value of parameter, interpreted as array of BYTE

      character
     & param*3,    ! Subroutine input argument (see above)
     & path*(*),   ! Subroutine input argument (see above)
     & cval*(*)    ! Subroutine input/output argument (see above)

      logical*4
     & comparing   ! Loop control during consistency check

      equivalence (i2val,bval)

      call get_opus_prm(luns,ipoint,param,iend,mrs,i2val,rs,typ)
      if(rs.lt.1) then
         errnum=4
         write(*,'(4a)')'Warning: ',
     &   path(1:lnbc(path)),' has no value for ',param
      elseif(typ.ne.exptyp) then
         errnum=4
         write(*,'(4a)')'Error: ',
     &   path(1:lnbc(path)),' has bad type for ',param
         if(verbose.ge.5) then
            write(*,*) param,'_type=',typ
         endif
      elseif(slicnt.eq.1) then
         do indexa=1,2*rs
            cval(indexa:indexa)=char(bval(indexa))
         enddo 
         if(verbose.ge.5) then
            write(*,*) param,'_raw=',(bval(indexb),indexb=1,2*rs)
         endif
         length=index(cval(1:2*rs),char(0))-1
c DG090212
c         if(length.eq.(-1)) length=2*rs  ! Work-around for missing terminator
         if(length.le.0) length=2*rs  ! Work-around for missing terminator or null string
c
         if(length.lt.minlen) then
            errnum=4
            write(*,'(5a)')'Error: ',
     &      path(1:lnbc(path)),' has incorrect ',param,' length.'
         elseif(verbose.ge.4) then
            write(*,*) param,'=',cval(1:length)
         endif
      else
         comparing=.true.
         indexa=1
         do while(comparing)
            if(bval(indexa).eq.0) then
               comparing=.false.
               length=indexa-1
            elseif((cval(indexa:indexa).ne.char(bval(indexa))).or.
     &      (indexa.eq.(2*rs))) then
               comparing=.false.
               errnum=4
               write(*,'(4a)')'Error: get_opus_string: ',
     &         path(1:lnbc(path)),' has inconsistent value for ',param
            else
               indexa=indexa+1
            endif
         enddo 
      endif

      return
      end
c ===================================================================
      subroutine get_opus_prm(luns,ipoint,param,iend,mrs,i2val,rs,typ)
c
c  NOTE: This soubroutine should not be used outside of this file.
c        The typed subroutines above should be used instead.
c
c  Reads an OPUS file and extracts the value of the parameter PARAM
c  from the parameter block starting at IPOINT.
c
c  Input:
c    luns      I*4    Logical Unit Number of file (already open)
c    ipoint    I*4    Pointer to parameter block (in byte units)
c    param     C*(3)  Name of parameter to search for
c    iend      I*4    Endian Flag
c    mrs       I*2    Maximum value for RS
c
c  Output:
c    i2val(mrs)I*2    Parameter value
c    rs        I*2    Length of parameter value (in 16-bit units)
c    typ       I*2    Type of the parameter, read from file
c
c On exit, rs<=0 indicates failure to find PARAM
c          rs>0  indicates success, with the length of data field = 2*RS bytes
c
      implicit none

      integer*4 
     & luns,       ! Subroutine input argument (see above)
     & ipoint,     ! Subroutine input argument (see above)
     & iend,       ! Subroutine input argument (see above)
     & bigendian,  ! Named constant for big endian detection 
     & reclen,     ! Record length in bytes for file reads
     & irec,       ! Record number during file reads
     & is,         ! Loop index during read of parameter value
     & i4var,      ! Local variable used for byte swapping
     & ierr,       ! Error status returned by file read
     & fnbc        ! Integer function First Non-Blank Character in string

      parameter (bigendian=1)
      parameter (reclen=2)

      integer*2 
     & mrs,        ! Subroutine input argument (see above)
     & i2val(mrs), ! Subroutine output argument (see above)
     & rs,         ! Subroutine output argument (see above)
     & typ,        ! Subroutine output argument (see above)
     & i2var(4)    ! Local variable used for byte swapping

      real*8
     & r8var       ! Local variable used for byte swapping

      character
     & param*3,    ! Subroutine input argument (see above)
     & head*4,     ! Current value of parameter name during search
     & stringa*11  ! String used to format integer display

      equivalence (i4var,i2var,r8var)
c
c  Search parameter block for PARAM label.
c
      irec=ipoint/reclen
      ierr=0
      rs=1
      head='....'
      do while ((head(:3).ne.param(:3)).and.(rs.gt.0)) 
         read(unit=luns,rec=1+irec,iostat=ierr) head(1:2)
         if(ierr.eq.0) then
            read(unit=luns,rec=2+irec,iostat=ierr) head(3:4)
         endif
         if(ierr.eq.0) then
            read(unit=luns,rec=3+irec,iostat=ierr) typ
         endif
         if(ierr.eq.0) then
            read(unit=luns,rec=4+irec,iostat=ierr) rs
         endif
         if(ierr.ne.0) then
            rs=0
            write(*,'(a)')
     &      'Error: file read failed during parameter search'
         elseif(iend.eq.bigendian) then
            call rbyte(rs,2,1)
            call rbyte(typ,2,1)
         endif
         irec=irec+4+rs
      enddo  ! while((head.ne.param).and.(rs.gt.0))
c
c  Check consistency of 'typ' and 'rs'.
c
      if(rs.le.0) then
         continue
      elseif((typ.eq.0).and.(rs.ne.2)) then
         write(stringa,'(i11)') rs
         write(*,'(4a)')'Error: INT32 parameter ',param,
     &   ' has bad length of ',stringa(fnbc(stringa):11)
         rs=0
      elseif((typ.eq.1).and.(rs.ne.4)) then
         write(stringa,'(i11)') rs
         write(*,'(4a)')'Error: REAL64 parameter ',param,
     &   ' has bad length of ',stringa(fnbc(stringa):11)
         rs=0
      elseif((typ.gt.4).or.(typ.lt.0)) then
         write(stringa,'(i11)') typ
         write(*,'(4a)')'Error: parameter ',param,
     &   ' has invalid type of ',stringa(fnbc(stringa):11)
         rs=0
      elseif(rs.gt.mrs) then
         write(stringa,'(i11)') rs
         write(*,'(2a)')'Error: increase parameter MRS to ',
     &   stringa(fnbc(stringa):11)
         rs=0
      endif
c
c  Read the value of PARAM.
c
      if(rs.gt.0) then
         irec=irec-rs
         do is=1,rs
            if(ierr.eq.0) then
               read(unit=luns,rec=irec+is,iostat=ierr) i2val(is)
            endif
            if(ierr.ne.0) then
               rs=0
               write(*,'(a)')
     &         'Error: file read failed during parameter read'
            endif
         enddo
c
c  Perform byte swapping for INT32 and REAL64 on big endian CPUs.
c
         if(ierr.ne.0) then
            continue
         elseif((typ.eq.0).and.(iend.eq.bigendian)) then
            do is=1,rs
               i2var(is)=i2val(is)
            enddo
            call rbyte(i4var,4,1)
            do is=1,rs
               i2val(is)=i2var(is)
            enddo
         elseif((typ.eq.1).and.(iend.eq.bigendian)) then
            do is=1,rs
               i2var(is)=i2val(is)
            enddo
            call rbyte(r8var,8,1)
            do is=1,rs
               i2val(is)=i2var(is)
            enddo
         endif
      endif     ! (rs.gt.0)
      return
      end
c ===================================================================
      subroutine test_opus_prm(luns,ipoint,param,iend,typ)
c
c  Reads an OPUS file and checks the type of the parameter PARAM
c  in the parameter block starting at IPOINT.
c
c  Input:
c    luns      I*4    Logical Unit Number of file (already open)
c    ipoint    I*4    Pointer to parameter block (in byte units)
c    param     C*(3)  Name of parameter to search for
c    iend      I*4    Endian Flag
c
c  Output:
c    typ       I*2    Type of the parameter, read from file
c
c On exit, typ=-1 indicates parameter not found or error
c          typ=0  indicates parameter is INT32 (I*4)
c          typ=1  indicates parameter is REAL64 (R*8)
c          typ=2  indicates parameter is STRING
c          typ=3  indicates parameter is ENUM
c          typ=4  indicates parameter is SENUM
c
      implicit none

      integer*4 
     & luns,       ! Subroutine input argument (see above)
     & ipoint,     ! Subroutine input argument (see above)
     & iend,       ! Subroutine input argument (see above)
     & bigendian,  ! Named constant for big endian detection 
     & reclen,     ! Record length in bytes for file reads
     & irec,       ! Record number during file reads
     & ierr,       ! Error status returned by file read
     & fnbc        ! Integer function First Non-Blank Character in string

      parameter (bigendian=1)
      parameter (reclen=2)

      integer*2 
     & typ,        ! Subroutine output argument (see above)
     & rs          ! Size of paramaters found in header block

      character
     & param*3,    ! Subroutine input argument (see above)
     & head*4,     ! Current value of parameter name during search
     & stringa*11  ! String used to format integer display
c
c  Search parameter block for PARAM label.
c
      irec=ipoint/reclen
      ierr=0
      rs=1
      head='....'
      do while ((head(:3).ne.param(:3)).and.(rs.gt.0)) 
         read(unit=luns,rec=1+irec,iostat=ierr) head(1:2)
         if(ierr.eq.0) then
            read(unit=luns,rec=2+irec,iostat=ierr) head(3:4)
         endif
         if(ierr.eq.0) then
            read(unit=luns,rec=3+irec,iostat=ierr) typ
         endif
         if(ierr.eq.0) then
            read(unit=luns,rec=4+irec,iostat=ierr) rs
         endif
         if(ierr.ne.0) then
            rs=0
            write(*,'(a)')
     &      'Error: file read failed during parameter search'
         elseif(iend.eq.bigendian) then
            call rbyte(rs,2,1)
            call rbyte(typ,2,1)
         endif
         irec=irec+4+rs
      enddo  ! while((head.ne.param).and.(rs.gt.0))
c
c  Check consistency of 'typ' and 'rs'.
c
      if(rs.le.0) then
         continue
      elseif((typ.eq.0).and.(rs.ne.2)) then
         write(stringa,'(i11)') rs
         write(*,'(2a)')'Error: INT32 parameter has bad length of ',
     &   stringa(fnbc(stringa):11)
         rs=0
      elseif((typ.eq.1).and.(rs.ne.4)) then
         write(stringa,'(i11)') rs
         write(*,'(2a)')'Error: REAL64 parameter has bad length of ',
     &   stringa(fnbc(stringa):11)
         rs=0
      elseif((typ.gt.4).or.(typ.lt.0)) then
         write(stringa,'(i11)') typ
         write(*,'(2a)')'Error: invalid parameter type of ',
     &   stringa(fnbc(stringa):11)
         rs=0
      endif
c
c  Return "not found" value if not found or if any error
c
      if(rs.le.0) then
         typ=1
         typ=-typ        ! Silly construction to avoid FTNCHEK warning
      endif
      return
      end
