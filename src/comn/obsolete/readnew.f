      subroutine readnew(ipath,mpath,id,nstart,ncount,ihead,ibuf,nip)
c
c  Subroutine to read MkIV interferogram files (binary, little endian)
c  by using the map file produced by 'maptelem' to locate the data
c  sections.  The map is an ASCII file; it starts with a pre-amble that
c  contains the date and time of the start of run, and the run direction.
c
c  The next section of the map consists of three columns:
c  1) The number of missing points before the current data segment.
c  2) The number of contiguous data points in the segment.
c  3) The number of bytes to skip in the file from the current read
c     position to get to the data segment.
c
c  Finally, a post-amble section gives the date and time of the end of
c  run, and repeats the run direction.
c
c  Each telemetry data frame contains 120 data values (60 HgCdTe + 60 InSb)
c  plus 8 words of housekeeping (sync words, frame tag, frame type, etc).
c
c  This subroutine returns the header and trailer information as straight
c  I*4 extracted from the map file.  The year includes the century.  The
c  run direction follows the simple convention that +1 means forward and
c  -1 means reverse.
c
c  There are two ways of using this subroutine:
c
c  1) If one sets 'nstart' to 1 and uses for 'ncount' the declared dimension
c     of array 'ibuf', the returned data will be the complete interferogram
c     and 'nip' will contain the length of the available data.  However
c     if the data buffer 'ibuf' is too small, the interferogram will be
c     truncated and 'nip' will be set to the size of the array.  This usage
c     of the subroutine is similar to that of the original version of
c     'readnew' and that of 'readigram'.
c
c  2) Any part of the interferogram may be read by setting 'nstart' and
c     'ncount' accordingly.  If the requested 'nstart' is past the end of
c     the interferogram, the subroutine will return no data and 'nip' will
c     be zero.  If the end of the request goes past the end of run, 'nip'
c     will be set to the number of points available hence truncating the
c     request.  A negative 'ncount' is illegal and will result in a stop
c     with an error message.  A count of zero is legal and will yield the
c     correct information in 'ihead' and will result in an 'nip' of zero.
c
c  NOTES: Unless there is a data gap, this subroutine will return the
c         same contents in 'ibuf' as the original 'readigram'.  Data
c         gaps take up the corresponding locations in 'ibuf', and these
c         are filled with alternating +/-2000 DN values to get the
c         attention of 'birdfix'.
c
c         The subroutine 'readigram' would warn about losses of sync,
c         but use the frame anyway.  The 'maptelem'+'readnew' combination
c         skips the frames with a loss of lock because their contents
c         cannot be trusted.  Hence the difference between the values of
c         'nip' returned by the two subroutines, is the size of the gaps
c         minus sixty times the number frames with loss of lock that had a
c         data frame id.
c
c  Input arguments:
c        IPATH      The path to the interferogram
c        MPATH      The path to the map file
c        ID         The detector of interest (HgCd=1; InSb=2)
c        NSTART     Index of the first interferogram point to be read
c        NCOUNT     Number of requested interferogram points
c
c  Output arguments:
c        IHEAD      The information from the telemetry header and trailer
c        IBUF       The array of interferogram values
c        NIP        The number of interferogram points read
c
c  Date: 10-Jan-2003,  JFB
c
      implicit none

      integer*4
     & id,             ! Subroutine input argument (see above)
     & nstart,         ! Subroutine input argument (see above)
     & ncount,         ! Subroutine input argument (see above)
     & ihead(16),      ! Subroutine output argument (see above)
     & nip,            ! Subroutine output argument (see above)
     & npos,           ! Position inside the full interferogram data
     & iendian,        ! Endianness of computer (0=big, 1=little)
     & krec,           ! Number of the frame to be read (record)
     & kcur,           ! Number of the frame currently in memory
     & indexa,         ! General loop index
     & ivalue,         ! General integer value
     & icount,         ! Number of points in data segment
     & igap,           ! Number of points in gap between data segments
     & iskip,          ! Number of bytes in non-data segment (e.g. filler)
     & frmpts,         ! Number of points per channel in a data frame
     & frmlen,         ! Number of I*2 data words per telemetry frame
     & frmbytes        ! Number of bytes per telemetry frame
      parameter (frmpts=60)
      parameter (frmlen=128)
      parameter (frmbytes=256)

      integer*2
     & ibuf(ncount),   ! Subroutine output argument (see above)
     & itmp(frmlen)    ! Buffer to hold the current data frame

      character
     & ipath*(*),      ! Subroutine input argument (see above)
     & mpath*(*),      ! Subroutine input argument (see above)
     & inpstr*40,      ! String used during map reading
     & rundir          ! Direction of run: F or R (forward or reverse)
c
c ---------------------------------------------------------------
c
c  This section of code performs all legality checks and opens the
c  two input files
c
c  Trap illegal input arguments
c
      if(ncount.lt.0) stop 'Readnew called with a negative count'
c
c  Open map file and check its format version number
c
      open(20,file=mpath,status='old',err=98)
      read(20,*,end=99)ivalue
      if(ivalue.ne.2) stop 'Map file has incorrect version number'
c
c  Check the instrument name and its data format
c
      read(20,*,end=99)inpstr,ivalue
      indexa=index(inpstr,'MkIV')
      if(indexa.ne.1) stop 'Map is for an unsupported instrument'
      if(ivalue.ne.1) stop 'Map is for an unsupported data format'
c
c  Open interferogram file
c
      open(19,file=ipath,status='old',access='direct',recl=frmbytes)
c
c ---------------------------------------------------------------
c
c  Now deal with the header information: read it from the pre-amble
c  of the map file
c
      read(20,*,end=99)ihead(1),ihead(2),ihead(3),ihead(4),
     &                 ihead(5),ihead(6),ihead(7),rundir
c
c  Convert the run direction information into an integer
c
      if(rundir.eq.'F') then
        ihead(8)=1
      elseif(rundir.eq.'R') then
        ihead(8)=-1
      else
        stop 'Invalid run direction in map file'
      endif
c
c ---------------------------------------------------------------
c
c  Skip the rest of the map file pre-amble (if any) and the delimiter
c
      read(20,*,end=99)inpstr
      do while (inpstr(1:1).ne.'=')
        read(20,*,end=99)inpstr
      end do
c
c ---------------------------------------------------------------
c
c  Next read the data from the interferogram file into the output buffer,
c  using the map file as a guide.
c
c  Initialize all loop variables
c
      nip=0      ! Number of points read so far
      npos=0     ! Positioned just before start of interferogram 
      krec=1     ! Start reading file at frame 1
      kcur=0     ! Force a disk read by positioning to non-existent frame
      indexa=id  ! Position frame content pointer to proper channel
      call endian(iendian)              !  SUN iendian=0;    PC iendian=1
c
c  Read next line of map file to initialize loop variable 'inpstr(1:1)'
c
      read(20,'(a)',end=99)inpstr
c
c  The loop will terminate either when the requested data length has been
c  met (test on 'nip') or when the end of the interferogram has been
c  reached (signaled by the end of the data section in the map file,
c  hence the test on 'inpstr(1:1)').
c 
      do while ((nip.lt.ncount).and.(inpstr(1:1).ne.'='))
c
c  Interpret one line of the map into the three integers
c
        read(inpstr,*)igap,icount,iskip
c
c  If there is a data gap inside the requested data, fill the corresponding
c  locations of the output buffer with a +/- 2000 DN oscillation.  Note that
c  0 DN has a value of 2048 in the telemetry format.  Negative gaps are
c  handled right after the while loop.
c
        do while ((igap.gt.0).and.(nip.lt.ncount))
          npos=npos+1
          if(npos.ge.nstart) then
            nip=nip+1
            ibuf(nip)=48+4000*mod(nip/60,2)
            if(iendian.eq.0) call byterev(ibuf(nip),2) ! In-place byte reversal
          endif
          igap=igap-1
        end do
        if(igap.lt.0) then
          npos=npos+igap
          if(npos.lt.0)stop'Tried to backup over start of interferogram' 
          nip=nip+igap
          if(nip.lt.0) nip=0 
        endif
c
c  Move the interferogram file-pointer pair (krec, indexa) by 'iskip' bytes
c
        krec=krec+(iskip/frmbytes)
        indexa=indexa+(mod(iskip,frmbytes)/2)
        if(indexa.gt.2*frmpts) then
          indexa=indexa-2*frmpts
          krec=krec+1
        elseif(indexa.lt.1) then
          indexa=indexa+2*frmpts
          krec=krec-1
        endif
        if(krec.lt.1) stop 'Tried to backup over start of file' 
c
c  Read the data points listed on this map line; this can span several frames
c
        do while ((icount.gt.0).and.(nip.lt.ncount))
c
c  See if we must read one telemetry data frame
c
          if((krec.ne.kcur).and.((npos+frmpts).ge.nstart)) then
            read(19,rec=krec)itmp
            kcur=krec
          endif
c
c  Transfer requested channel from 'itmp' to 'ibuf'
c
          do while((indexa.le.(2*frmpts)).and.(icount.gt.0).and.
     &             (nip.lt.ncount))
            npos=npos+1
            if(npos.ge.nstart) then
              nip=nip+1
              ibuf(nip)=itmp(indexa)
            endif
            indexa=indexa+2
            icount=icount-1
          end do
          if(indexa.gt.2*frmpts) then
            indexa=indexa-2*frmpts
            krec=krec+1
          endif
c
        end do        ! while ((icount.gt.0).and.(nip.lt.ncount))
c
c  Read the next map line so that it is ready for the while test
c
        read(20,'(a)',end=99)inpstr
c
      end do          ! while ((nip.lt.ncount).and.(inpstr(1:1).ne.'='))
c
c ---------------------------------------------------------------
c
c  Now deal with the trailer information.  But first make sure that
c  we are no longer in the data section:
c
      do while (inpstr(1:1).ne.'=')
        read(20,*,end=99)inpstr
      end do
c
c  Read trailer info from the post-amble of the map file
c
      read(20,*,end=99)ihead(9),ihead(10),ihead(11),ihead(12),
     &                 ihead(13),ihead(14),ihead(15),rundir
c
c  Convert the run direction information into an integer
c
      if(rundir.eq.'F') then
        ihead(16)=1
      elseif(rundir.eq.'R') then
        ihead(16)=-1
      else
        stop 'Invalid run direction in map file'
      endif
c
c ---------------------------------------------------------------
c
c  Finally, we deal with then endianness issue: because the MkIV telemetry
c  file is little-endian, big-endian computers must byte swap when reading 
c  these data.
c
      if(iendian.eq.0) then
        do indexa=1,nip
          call byterev(ibuf(indexa),2)  ! In-place byte reversal
        end do
      endif
c
c ---------------------------------------------------------------
c
c  All done !  Close the two input files
c
      close(19)
      close(20)
      return
c
98    stop 'Cannot open map file'
99    stop 'Unexpected end of map file'
      end
