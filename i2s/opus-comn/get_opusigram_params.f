      subroutine get_opusigram_params(path,verbose,msl,mch,mi4,mr8,
     & chan1,chan2,errnum,slicnt,inrun,nptvec,bpdata,timsli,nss,
     & tpx,i4head,r8head,DTCstr,INSstr)
c
c  Input:
c    path      C*(*)  Path to IFS125 interferogram slice file
c    verbose   I*4    Level of verbosity for displayed messages
c    msl       I*4    Maximum number of interferogram slices per scan set
c    mch       I*4    Maximum number of data channels
c    mi4       I*4    Maximum number of I*4 items in file header
c    mr8       I*4    Maximum number of R*8 items in file header
c    chan1     I*4    Starting channel number to process
c    chan2     I*4    Ending channel number to process
c
c  Input/output:
c    errnum    I*4    Error code (0=ok, <0=fatal, >0=recoverable)
c    slicnt    I*4    Slice counter, points into nptvec, bpdata, and timsli
c
c  Output:
c    inrun     L*4    Set to false at end of NSS runs
c    nptvec(msl)I*4   Number of PoinTs, from OPUS header
c    bpdata(msl,mch)I*4 Byte pointers into the data blocks of slices
c    timsli(msl)R*8   Time vector contains one entry for each slice
c    nss       I*4    Number of Sample Scans in this set of slices
c    ass       I*4    Number of Sample Scans in this set of slices
c    tpx       I*4    Number of points in one FWD+REV pair
c    i4head(mi4)I*4   Vector to hold the I*4 header items
c    r8head(mr8)R*8   Vector to hold the R*8 header items
c    DTCstr     C*40  Variable to hold detector description
c    INSstr     C*40  Variable to hold instrument description
c
      implicit none

      integer*4 ii,
     & verbose,    ! Subroutine input argument (see above)
     & msl,        ! Subroutine input argument (see above)
     & mch,        ! Subroutine input argument (see above)
     & mi4,        ! Subroutine input argument (see above)
     & mr8,        ! Subroutine input argument (see above)
     & chan1,      ! Subroutine input argument (see above)
     & chan2,      ! Subroutine input argument (see above)
     & errnum,     ! Subroutine input/output argument (see above)
     & slicnt,     ! Subroutine input/output argument (see above)
     & nptvec(msl),! Subroutine output argument (see above)
     & bpdata(msl,mch),! Subroutine output argument (see above)
     & nss,        ! Subroutine output argument (see above)
     & ass,        ! Subroutine output argument (see above)
     & tpx,        ! Subroutine output argument (see above)
     & i4head(mi4),! Subroutine output argument (see above)
     & iendian,    ! Endianness of computer
     & bigendian,  ! Named constant for big endian detection 
     & idtagd(2),  ! Block type for the two channel data, master/slave
     & idtags(2),  ! Block type for the two channel data status, master/slave
     & npts(2),    ! Number of points in this slice for the two detectors
     & inpoint(2), ! Pointers to data values in slice for the two detectors
     & ipdstat(2), ! Pointers to data status parameter blocks, slave/master
     & ipinstr,    ! Pointer to instrument status parameter block
     & ipaqpar,    ! Pointer to standard acquisition parameter block
     & ipoptpar,   ! Pointer to optics parameter block
     & iporgpar,   ! Pointer to sample origin parameter block
     & channel,    ! Channel number (1=slave, 2=master)
     & inpstat,    ! Value of IOSTAT from all file open/read/close
     & dstatcnt,   ! Count of data status blocks that were found
     & luns,       ! Logical Unit Number for slice input file
     & reclen,     ! Record length in bytes for single values
     & magicin,    ! Value of first word from the slice file
     & ifilver(2), ! Integer EQUIVALENT to double-precision 'filverin'
     & ipdb,       ! Byte pointer into directory block
     & mndb,       ! Maximum number of directory entries
     & ndb,        ! Number of directory entries
     & indexa,     ! General loop index
     & blocktyp,   ! Directory entry for the type of data/parameter
     & blocklen,   ! Directory entry for the length of data/parameter
     & blockpnt,   ! Directory entry for the pointer to data/parameter
     & fnbc,       ! Integer function First Non-Blank Character in string
     & lnbc,       ! Integer function Last Non-Blank Character in string
     & strlen,     ! String length returned from parameter read
     & timlen,     ! String length returned from time read
     & datlen,     ! String length returned from date read
     & iy,         ! Four-digit year
     & im,         ! Month
     & id,         ! Day
     & hh,         ! Hour
     & mm,         ! Minute
     & ss,         ! Second
     & ms,         ! Millisecond
     & jd,         ! Julian Day
     & j1,         ! Julian Day of 1-Jan-2000
     & fxv,        ! First X Value, from OPUS header
     & lxv,        ! Last X Value, from OPUS header
     & fxvchk,     ! Used only during consistency check of FXV
     & lxvchk      ! Used only during consistency check of LXV

      parameter (bigendian = 1)

      integer*2
     & paramtyp    ! Type of OPUS parameters, used to check for existence

      real*8
     & timsli(msl),! Subroutine output argument (see above)
     & r8head(mr8),! Subroutine output argument (see above)
     & filverin,   ! Value of file version read in
     & r8fxv,      ! Temporary variable used in reading FXV
     & r8lxv,      ! Temporary variable used in reading LXV
     & velocity(10)! List of velocity values used with older IFS120
ccc     & datarate    ! Interferogram sampling rate in Hz

      character
     & path*(*),   ! Subroutine input argument (see above)
     & cval*99,    ! Value of string parameter
     & timval*99,  ! Value of time string
     & datval*99,  ! Value of date string
     & chan(2)*6,  ! Identifiers for each channel (detector name)
     & stringa*11, ! String used to format integer display
     & stringb*11, ! String used to format integer display
     & ch1*1,      ! Character used in format checks
     & ch2*1,      ! Character used in format checks
     & ch3*1,      ! Character used in format checks
     & DTCstr*(*), ! Detector description string
     & INSstr*(*)  ! Instrument description string

      logical*4
     & inrun       ! Subroutine output argument (see above)

      include 'opus_constants.inc'
      include 'header_indices.inc'

      parameter (luns=21)

      equivalence (ifilver,filverin)

      data chan /'slave ','master'/
      data velocity/0.d0,0.d0,4.d0,5.d0,7.5d0,10.d0,20.d0,40.d0,
     & 60.d0,80.d0/
c
c  Initialize variables.
c
      call getendian(iendian)

      idtagd(2)=dbbampl+dbbsamp+dbbigrm  ! Master interferogram data
      idtagd(1)=idtagd(2)+dbbslav        ! Slave interferogram data
      idtags(2)=idtagd(2)+dbbdstat       ! Master data status
      idtags(1)=idtagd(1)+dbbdstat       ! Slave data status
c
c  Indicate that parameter blocks have not yet been found.
c
      do channel=chan1,chan2
        inpoint(channel)=0
        ipdstat(channel)=0
      enddo
      dstatcnt=0
      ipinstr=0
      ipaqpar=0
      ipoptpar=0
      iporgpar=0
c
c  This subroutine is divided into two major sections: the slice file
c  is first opened with a record length of 4 to search the directory
c  block and locate the data and parameter blocks. 
c
      if(errnum.eq.0) then
        reclen=4              ! Four bytes to read directory blocks
        open(unit=luns,file=path,form='unformatted',status='old',
     &   access='direct',recl=reclen,iostat=inpstat) 
        if(inpstat.ne.0) then
          errnum=4
          write(*,'(2a)')
     &     'Error: open failed on file ',path(1:lnbc(path))
        endif
c
c  Read and check OPUS magic 4-byte integer.
c
        if(errnum.eq.0) then
          read(unit=luns,rec=1,iostat=inpstat)magicin
          if(inpstat.ne.0) then
            errnum=4
            write(*,'(2a)')
     &       'Error: magic read failed on file ',path(1:lnbc(path))
          else
            if(iendian.eq.bigendian) call i4rev(magicin)
            if(verbose.ge.5) write(*,*) 'File_magic=',magicin
            if(magicin.ne.magic) then
              errnum=4
              write(*,'(3a)')
     &         'Error: ',path(1:lnbc(path)),' has bad magic'
            endif
          endif
        endif
c
c  Read and check OPUS version 8-byte floating-point.
c
        do indexa=1,2
          if(errnum.eq.0) then
            read(unit=luns,rec=1+indexa,iostat=inpstat)ifilver(indexa)
            if(inpstat.ne.0) then
              errnum=4
              write(*,'(2a)')
     &         'Error: version read failed on file ',path(1:lnbc(path))
            endif
          endif
        enddo
        if(errnum.eq.0) then
          if(iendian.eq.bigendian) call r8rev(filverin)
          if(verbose.ge.5) write(*,*) 'File_version=',filverin
          if(dabs(filverin-filevers).gt.0.01d0) then
            errnum=4
            write(*,'(3a)')
     &       'Error: ',path(1:lnbc(path)),' has bad version'
          endif
        endif
c
c  Read and check byte pointer into directory block.
c
        if(errnum.eq.0) then
          read(unit=luns,rec=4,iostat=inpstat)ipdb
          if(inpstat.ne.0) then
            errnum=4
            write(*,'(2a)')
     &       'Error: read failed at directory pointer, file ',
     &       path(1:lnbc(path))
          else
            if(iendian.eq.bigendian) call i4rev(ipdb)
            if(verbose.ge.5) write(*,*) 'Directory_pointer=',ipdb
            if(ipdb.lt.0) then
              errnum=4
              write(*,'(2a)')
     &         'Error: directory pointer is negative, file ',
     &         path(1:lnbc(path))
            elseif(ipdb.lt.24) then
              errnum=4
              write(*,'(2a)')
     &         'Error: directory block overlaps header, file ',
     &         path(1:lnbc(path))
            elseif(mod(ipdb,reclen).ne.0) then
              errnum=4
              write(*,'(2a)')
     &         'Error: directory block not on a word boundary, file ',
     &         path(1:lnbc(path))
            endif
          endif
        endif
c
c  Read and check maximum size of directory.
c
        if(errnum.eq.0) then
          read(unit=luns,rec=5,iostat=inpstat)mndb
          if(inpstat.ne.0) then
            errnum=4
            write(*,'(2a)')
     &       'Error: read failed at max directory size, file ',
     &       path(1:lnbc(path))
          else
            if(iendian.eq.bigendian) call i4rev(mndb)
            if(verbose.ge.5) write(*,*) 'Directory_max_size=',mndb
            if(mndb.lt.0) then
              errnum=4
              write(*,'(2a)')
     &         'Error: directory max size is negative, file ',
     &         path(1:lnbc(path))
            endif
          endif
        endif
c
c  Read and check current size of directory (i.e., number of entries).
c
        if(errnum.eq.0) then
          read(unit=luns,rec=6,iostat=inpstat)ndb
          if(inpstat.ne.0) then
            errnum=4
            write(*,'(2a)')
     &       'Error: read failed at directory size, file ',
     &       path(1:lnbc(path))
          else
            if(iendian.eq.bigendian) call i4rev(ndb)
            if(verbose.ge.5) write(*,*) 'Directory_size=',ndb
            if(ndb.lt.0) then
              errnum=4
              write(*,'(2a)')
     &         'Error: directory size is negative, file ',
     &         path(1:lnbc(path))
            elseif(ndb.gt.mndb) then
              errnum=4
              write(*,'(2a)')
     &         'Error: directory size exceeds maximum, file ',
     &         path(1:lnbc(path))
            endif
          endif
        endif
c
c  Convert pointer to directory block from bytes to records.
c
        if(errnum.eq.0) then
          ipdb=ipdb/reclen
c
c  Search directory entries for channel data blocks,
c  channel data status blocks, and other parameter blocks.
c
          do indexa=1,3*ndb,3

            if(errnum.eq.0) then
              read(unit=luns,rec=indexa+ipdb,iostat=inpstat) blocktyp
              if(inpstat.ne.0) then
                errnum=4
                write(*,'(2a)')
     &           'Error: read failed at block type, file ',
     &           path(1:lnbc(path))
              else
                if(iendian.eq.bigendian) call i4rev(blocktyp)
                blocktyp=mod(blocktyp,2**30)   ! Some old OPUS files need this
                if(verbose.ge.5) write(*,*) 'Block_type=',blocktyp
              endif
            endif

            if(errnum.eq.0) then
              read(unit=luns,rec=indexa+ipdb+1,iostat=inpstat) blocklen
              if(inpstat.ne.0) then
                errnum=4
                write(*,'(2a)')
     &           'Error: read failed at block length, file ',
     &           path(1:lnbc(path))
              else
                if(iendian.eq.bigendian) call i4rev(blocklen)
                if(verbose.ge.5) write(*,*) 'Block_length=',blocklen
              endif
            endif

            if(errnum.eq.0) then
              read(unit=luns,rec=indexa+ipdb+2,iostat=inpstat) blockpnt
              if(inpstat.ne.0) then
                errnum=4
                write(*,'(2a)')
     &           'Error: read failed at block pointer, file ',
     &           path(1:lnbc(path))
              else
                if(iendian.eq.bigendian) call i4rev(blockpnt)
                if(verbose.ge.5) write(*,*) 'Block_pointer=',blockpnt
              endif
            endif

            if(errnum.eq.0) then
              do channel=chan1,chan2
                if(blocktyp.eq.idtagd(channel)) then
                  npts(channel)=blocklen
                  inpoint(channel)=blockpnt
                  if(mod(inpoint(channel),reclen).ne.0)then
                    errnum=4
                    write(*,'(3a)')'Error: ',
     &               chan(channel)(1:lnbc(chan(channel))),
     &               ' data block not on a word boundary'
                  endif
                elseif(blocktyp.eq.idtags(channel)) then
                  ipdstat(channel)=blockpnt
                  if(mod(ipdstat(channel),reclen).ne.0)then
                    errnum=4
                    write(*,'(4a)')'Error: ',
     &               chan(channel)(1:lnbc(chan(channel))),
     &               ' data status ',
     &               'parameter block not on a word boundary'
                  endif
                endif
              enddo
              if(blocktyp.eq.dbbaqpar) then
                ipaqpar=blockpnt
                if(mod(ipaqpar,reclen).ne.0)then
                  errnum=4
                  write(*,'(2a)')'Error: acquisition ',
     &             'parameter block not on a word boundary'
                endif
              elseif(blocktyp.eq.dbbinstr) then
                inrun=.false.
                ipinstr=blockpnt
                if(mod(ipinstr,reclen).ne.0)then
                  errnum=4
                  write(*,'(2a)')'Error: instrument status ',
     &             'parameter block not on a word boundary'
                endif
              elseif(blocktyp.eq.dbboptpar) then
                inrun=.false.
                ipoptpar=blockpnt
                if(mod(ipoptpar,reclen).ne.0)then
                  errnum=4
                  write(*,'(2a)')'Error: optics ',
     &             'parameter block not on a word boundary'
                endif
              elseif(blocktyp.eq.dbborgpar) then
                iporgpar=blockpnt
                if(mod(iporgpar,reclen).ne.0)then
                  errnum=4
                  write(*,'(2a)')'Error: sample origin ',
     &             'parameter block not on a word boundary'
                endif
              endif
            endif
          enddo  ! indexa=1,3*ndb,3
        endif

        close(unit=luns,iostat=inpstat)
        if(inpstat.ne.0) then
          errnum=4
          write(*,'(2a)')'Error: close failed on file ',path
        endif
      endif
c
c  Check that data and parameter blocks were found.  All slices
c  but the last one should contain a data block, a data status
c  parameter block and an acquisition parameter block.  The data
c  status parameters include NPT and DAT+TIM, so we need to increment
c  'slicnt' when that block is found.  'Slicnt' starts at 0 and is
c  used as a pointer into the NPT and time vectors, so the increment
c  must be done before saving NPT and time.  'Slicnt' also points
c  into the 'bpdata' (byte pointer to data) vector.
c
c  The very last slice contains the instrument status, optics, and
c  sample origin parameters, and it has no data, data status, or
c  acquisition parameters.  Since this last slice sets 'inrun' to
c  false, we can use that to guard the checks.
c
c  Special care is needed for the case of the single slice scan,
c  which contains data and all six parameter blocks: 'inrun' is
c  already set to false, and we must force the checks if the slice
c  counter 'slicnt' is at 0.
c
      if((errnum.eq.0).and.(inrun.or.(slicnt.eq.0))) then
        do channel=chan1,chan2
          if(inpoint(channel).eq.0) then
            errnum=4
            write(*,'(5a)')'Error: ',path(1:lnbc(path)),
     &       ' has no ',chan(channel)(1:lnbc(chan(channel))),
     &       ' data'
          endif
          if(ipdstat(channel).eq.0) then
            errnum=4
            write(*,'(5a)')'Error: ',path(1:lnbc(path)),
     &       ' has no ',chan(channel)(1:lnbc(chan(channel))),
     &       ' data status parameters'
          endif
        enddo
        if((errnum.eq.0).and.(slicnt.eq.msl)) then
          errnum=4
          write(*,'(a)')'Error: increase parameter MSL'
        else
          slicnt=slicnt+1
        endif
        if(ipaqpar.eq.0) then
          errnum=4
          write(*,'(3a)')'Error: ',
     &     path(1:lnbc(path)),' has no acquisition parameters'
        endif
      endif
      if((errnum.eq.0).and.(.not.inrun)) then
        if(ipinstr.eq.0) then
          errnum=4
          write(*,'(3a)')'Error: ',
     &     path(1:lnbc(path)),' has no instrument status parameters'
        endif
        if(ipoptpar.eq.0) then
          errnum=4
          write(*,'(3a)')'Error: ',
     &     path(1:lnbc(path)),' has no optics parameters'
        endif
        if(iporgpar.eq.0) then
          errnum=4
          write(*,'(3a)')'Error: ',
     &     path(1:lnbc(path)),' has no sample origin parameters'
        endif
      endif
c
c  This is the second major section of the subroutine: now that the
c  parameter blocks have been found, the slice file is re-opened with
c  a record length of 2 to obtain values of all parameters.
c
      if(errnum.eq.0) then
        reclen=2              ! Two bytes to read parameter blocks
        open(unit=luns,file=path,form='unformatted',status='old',
     &   access='direct',recl=reclen,iostat=inpstat) 
        if(inpstat.ne.0) then
          errnum=4
          write(*,'(2a)')'Error: open failed on file ',path
        endif
c
c  Read from acquisition parameter block.
c
        if((errnum.eq.0).and.(ipaqpar.ne.0)) then
          call get_opus_i4(luns,ipaqpar,'NSS',iendian,path,verbose,
     &     slicnt,errnum,nss)

          cval=char(0)
          if(i4head(i_aqmcode).eq.aqm_sd) then
            write(cval,'(a)') 'SD'
          elseif(i4head(i_aqmcode).eq.aqm_sf) then
            write(cval,'(a)') 'SF'
          else
            write(cval,'(a)') '  '
          endif
          call get_opus_string(luns,ipaqpar,'AQM',typ_enum,2,iendian,
     &     path,verbose,slicnt,errnum,strlen,cval)
          if(errnum.eq.0) then
            i4head(i_aqmcode)=0
            if(index(cval(1:strlen),'SD').eq.1)
     &       i4head(i_aqmcode)=aqm_sd
            if(index(cval(1:strlen),'SF').eq.1)
     &       i4head(i_aqmcode)=aqm_sf
            if(index(cval(1:strlen),'SN').eq.1)
     &       i4head(i_aqmcode)=aqm_sn
            if(i4head(i_aqmcode).eq.0) then
              errnum=4
              write(*,'(4a)')'Error: ',path(1:lnbc(path)),
     &         ' has unexpected AQM of ',cval(1:strlen)
            elseif(verbose.ge.5) then
              write(*,*) 'AQM_cooked=',i4head(i_aqmcode)
            endif
          endif   !  if((errnum.eq.0).and.(ipaqpar.ne.0))

          call test_opus_prm(luns,ipaqpar,'SGN',iendian,paramtyp)
          if(paramtyp.eq.typ_enum) then
            cval=char(0)
            write(stringa,'(i11)') i4head(i_sgnb)
            write(cval,'(a)') stringa(fnbc(stringa):11)
            call get_opus_string(luns,ipaqpar,'SGN',typ_enum,1,iendian,
     &       path,verbose,slicnt,errnum,strlen,cval)
            if(errnum.eq.0) then
              read(cval(1:strlen),*,iostat=inpstat) i4head(i_sgnb)
              if(inpstat.ne.0) then
                errnum=4
                write(*,'(3a)')'Error: ',
     &           path(1:lnbc(path)),' has incorrect SGN format'
              elseif(verbose.ge.5) then
                write(*,*) 'SGN_cooked=',i4head(i_sgnb)
              endif
            endif
          else
            i4head(i_sgnb)=(-1)
          endif

          call test_opus_prm(luns,ipaqpar,'SG2',iendian,paramtyp)
          if(paramtyp.eq.typ_enum) then
            cval=char(0)
            write(stringa,'(i11)') i4head(i_sgna)
            write(cval,'(a)') stringa(fnbc(stringa):11)
            call get_opus_string(luns,ipaqpar,'SG2',typ_enum,1,iendian,
     &       path,verbose,slicnt,errnum,strlen,cval)
            if(errnum.eq.0) then
              read(cval(1:strlen),*,iostat=inpstat) i4head(i_sgna)
              if(inpstat.ne.0) then
                errnum=4
                write(*,'(3a)')'Error: ',
     &           path(1:lnbc(path)),' has incorrect SG2 format'
              elseif(verbose.ge.5) then
                write(*,*) 'SG2_cooked=',i4head(i_sgna)
              endif
            endif
          else
            i4head(i_sgna)=(-1)
          endif

        endif        ! if((errnum.eq.0).and.(ipaqpar.ne.0))
c
c  Read from data status parameter blocks.  Because there is one such
c  block per channel, we use a do loop to scan all the requested
c  channels.  Most parameters (except MNY and MXY) will match between
c  channels.  This is checked internally by the "get_opus_XXX" routines,
c  provided that the argument after 'verbose' is not equal to 1.  Hence
c  the use of '1+dstatcnt' or 'slicnt+dstatcnt' below.  The latter form
c  is used for parameters that should also match between slices in
c  addition to between channels.
c
        do channel=chan1,chan2
          if((errnum.eq.0).and.(ipdstat(channel).ne.0)) then

            call get_opus_i4(luns,ipdstat(channel),'NPT',iendian,path,
     &       verbose,1+dstatcnt,errnum,nptvec(slicnt))

            call test_opus_prm(luns,ipdstat(channel),'TPX',iendian,
     &       paramtyp)
            if(paramtyp.eq.typ_i4) then
              call get_opus_i4(luns,ipdstat(channel),'TPX',iendian,path,
     &         verbose,slicnt+dstatcnt,errnum,tpx)
            else
              tpx=nptvec(slicnt)
            endif

            call get_opus_r8(luns,ipdstat(channel),'FXV',iendian,path,
     &       verbose,1+dstatcnt,errnum,r8fxv)
            fxv=int(r8fxv)

            call get_opus_r8(luns,ipdstat(channel),'LXV',iendian,path,
     &       verbose,1+dstatcnt,errnum,r8lxv)
            lxv=int(r8lxv)

            call get_opus_string(luns,ipdstat(channel),'DAT',typ_string,
     &       10,iendian,path,verbose,1+dstatcnt,errnum,datlen,datval)

            call get_opus_string(luns,ipdstat(channel),'TIM',typ_string,
     &       8,iendian,path,verbose,1+dstatcnt,errnum,timlen,timval)

            dstatcnt=dstatcnt+1
          endif
        enddo

        if((errnum.eq.0).and.(dstatcnt.gt.0)) then
          if(tpx.gt.0) then
            continue
          elseif((tpx.eq.(-1)).and.(slicnt.eq.1).and.(nss.eq.1)) then
             tpx=nptvec(slicnt)
             if(verbose.ge.3) then
              write(*,'(2a)')path(1:lnbc(path)),
     &         ' is a single-slice scan, using NPT for TPX'
             endif
          else
             errnum=4
             write(*,'(3a)')'Error: ',
     &       path(1:lnbc(path)),' has bad TPX value'
c   Vanessa's "TERRIBLE KLUDGE" to allow processing of lamp GPR spectra
             write(*,*)' Assuming single-slice scan. Using NPT for TPX'
             tpx=nptvec(slicnt)
             errnum=0
          endif
        endif

        if((errnum.eq.0).and.(dstatcnt.gt.0)) then
          if(datval(5:5).eq.'/') then      !  YYYY/MM/DD
            read(datval,'(i4,a1,i2,a1,i2)',iostat=inpstat)
     &       iy,ch1,im,ch2,id
          else                           !  DD/MM/YYYY
            read(datval,'(i2,a1,i2,a1,i4)',iostat=inpstat)
     &       id,ch1,im,ch2,iy
          endif
          if((inpstat.ne.0).or.(ch1.ne.'/').or.(ch2.ne.'/')) then
            errnum=4
            write(*,'(3a)')'Error: ',
     &       path(1:lnbc(path)),' has incorrect DAT format'
          elseif((iy.lt.2000).or.(im.lt.1).or.(im.gt.12)
     &                       .or.(id.lt.1).or.(id.gt.31)) then
            errnum=4
            write(*,'(3a)')'Error: ',
     &       path(1:lnbc(path)),' has incorrect DAT value'
          elseif(verbose.ge.5) then
            write(*,*) 'DAT_cooked=',iy,im,id
          endif
        endif
        if((errnum.eq.0).and.(dstatcnt.gt.0)) then
          call julian(iy,im,id,jd)
          call julian(2000,1,1,j1)
          jd=jd-j1
        endif

        if((errnum.eq.0).and.(dstatcnt.gt.0)) then
          if((timlen.gt.8).and.(timval(9:9).eq.'.')) then
            read(timval,'(i2,a1,i2,a1,i2,a1,i3)',iostat=inpstat)
     &       hh,ch1,mm,ch2,ss,ch3,ms
          else
            read(timval,'(i2,a1,i2,a1,i2)',iostat=inpstat)
     &       hh,ch1,mm,ch2,ss
            ms=0
            ch3='.'
          endif
          if((inpstat.ne.0).or.
     &     (ch1.ne.':').or.(ch2.ne.':').or.(ch3.ne.'.')) then
            errnum=4
            write(*,'(3a)')'Error: ',
     &       path(1:lnbc(path)),' has incorrect TIM format'
          elseif((hh.lt.0).or.(hh.gt.23).or.
     &           (mm.lt.0).or.(mm.gt.59).or.
     &           (ss.lt.0).or.(ss.gt.59).or.
     &           (ms.lt.0).or.(ms.gt.999)) then
            errnum=4
            write(*,'(3a)')'Error: ',
     &       path(1:lnbc(path)),' has incorrect TIM value'
          elseif(verbose.ge.5) then
            write(*,*) 'TIM_cooked=',hh,mm,ss,ms
          endif
        endif
        if((errnum.eq.0).and.(dstatcnt.gt.0)) then
          timsli(slicnt)=dble(((jd*24+hh)*60+mm)*60+ss)+
     &     (dble(ms)/1000.d0)
        elseif(dstatcnt.gt.0) then
          timsli(slicnt)=0.d0
        endif
c
c  Consistency checks:
c    1) The count of points in the data blocks should match NPT from header.
c    2) The first and last X values should match the cumulated NPT vector.
c  We only display one error message to avoid clutter.
c
        if((errnum.eq.0).and.(dstatcnt.gt.0)) then
          do channel=chan1,chan2
            if(inpoint(channel).gt.0) then
              if(npts(channel).ne.nptvec(slicnt)) then
                errnum=4
                write(stringa,'(i11)') nptvec(slicnt)
                write(stringb,'(i11)') npts(channel)
                write(*,'(8a)')'Error: ',path(1:lnbc(path)),
     &           ' has mismatch between NPT=',stringa(fnbc(stringa):11),
     &           ' and ',chan(channel)(1:lnbc(chan(channel))),
     &           ' point count of ',stringb(fnbc(stringb):11)
              else
                bpdata(slicnt,channel)=inpoint(channel)
              endif
            endif
          enddo
          if(errnum.eq.0) then
            fxvchk=0
            do indexa=1,slicnt-1
              fxvchk=fxvchk+nptvec(indexa)
            enddo
            fxvchk=mod(fxvchk,tpx)
            lxvchk=fxvchk+nptvec(slicnt)-1
            if(fxv.ne.fxvchk) then
              errnum=4
              write(stringa,'(i11)') fxv
              write(stringb,'(i11)') fxvchk
              write(*,'(6a)')'Error: ',path(1:lnbc(path)),
     &         ' has mismatch between FXV=',stringa(fnbc(stringa):11),
     &         ' and starting count of ',stringb(fnbc(stringb):11)
            elseif(lxv.ne.lxvchk) then
              errnum=4
              write(stringa,'(i11)') lxv
              write(stringb,'(i11)') lxvchk
              write(*,'(6a)')'Error: ',path(1:lnbc(path)),
     &         ' has mismatch between LXV=',stringa(fnbc(stringa):11),
     &         ' and ending count of ',stringb(fnbc(stringb):11)
            endif
          endif
        endif
c
c  Read from optics parameter block.
c
        if((errnum.eq.0).and.(ipoptpar.ne.0)) then
          cval=char(0)
          call get_opus_string(luns,ipoptpar,'VEL',typ_enum,1,iendian,
     &     path,verbose,1,errnum,strlen,cval)
          if(errnum.eq.0) then
            indexa=index(cval(1:strlen),'.')
            if(indexa.gt.0) then  ! new IFS12x: read velocity directly
              read(cval(1:strlen),*,iostat=inpstat)r8head(i_laserate)
              if(inpstat.ne.0) then
                errnum=4
                write(*,'(3a)')'Error: ',
     &           path(1:lnbc(path)),' has incorrect VEL format'
              endif
            else                  ! old IFS120: read velocity index
              read(cval(1:strlen),*,iostat=inpstat)indexa
              if(inpstat.ne.0) then
                errnum=4
                write(*,'(3a)')'Error: ',
     &           path(1:lnbc(path)),' has incorrect VEL index format'
              elseif((indexa.lt.1).or.(indexa.gt.10)) then
                errnum=4
                write(*,'(3a)')'Error: ',
     &           path(1:lnbc(path)),' has VEL index out of bounds'
              else
                r8head(i_laserate)=velocity(indexa)
              endif
            endif
          endif
          if(errnum.eq.0) then
            if(verbose.ge.5) then
              write(*,*) 'VEL_cooked=',r8head(i_laserate)
            endif
            r8head(i_laserate)=r8head(i_laserate)*1000.d0
          endif

          cval=char(0)
          call get_opus_string(luns,ipoptpar,'SRC',typ_senum,3,iendian,
     &     path,verbose,1,errnum,strlen,cval)
          if(errnum.eq.0) then
            i4head(i_srccode)=0
            if(index(cval(1:strlen),'Emission').eq.1)
     &       i4head(i_srccode)=src_sun
            if(index(cval(1:strlen),'Solar').eq.1)
     &       i4head(i_srccode)=src_sun
            if(index(cval(1:strlen),'Sun').eq.1)
     &       i4head(i_srccode)=src_sun
            if(index(cval(1:strlen),'NIR').gt.0)
     &       i4head(i_srccode)=src_nir
            if(index(cval(1:strlen),'MIR').gt.0)
     &       i4head(i_srccode)=src_mir
            if(index(cval(1:strlen),'Off All').eq.1)
     &       i4head(i_srccode)=src_off
            if(i4head(i_srccode).eq.0) then
              errnum=4
              write(*,'(4a)')'Error: ',path(1:lnbc(path)),
     &         ' has unexpected SRC of ',cval(1:strlen)
            elseif(verbose.ge.5) then
              write(*,*) 'SRC_cooked=',i4head(i_srccode)
            endif
          endif

          cval=char(0)
          call get_opus_string(luns,ipoptpar,'APT',typ_senum,4,iendian,
     &     path,verbose,1,errnum,strlen,cval)
          if(errnum.eq.0) then
            read(cval(1:index(cval,'mm')-1),*,iostat=inpstat)
     &       r8head(i_aptval)
            if(inpstat.ne.0) then
              errnum=4
              write(*,'(3a)')'Error: ',
     &         path(1:lnbc(path)),' has incorrect APT format'
            elseif(verbose.ge.5) then
              write(*,*) 'APT_cooked=',r8head(i_aptval)
            endif
          endif

          call test_opus_prm(luns,ipoptpar,'BMS',iendian,paramtyp)
          if(paramtyp.eq.typ_senum) then
            call get_opus_string(luns,ipoptpar,'BMS',typ_senum,3,
     &       iendian,path,verbose,1,errnum,strlen,cval)
            if(errnum.eq.0) then
              i4head(i_bmscode)=0
              if(index(cval(1:strlen),'CaF2').eq.1)
     &         i4head(i_bmscode)=bms_caf2
              if(index(cval(1:strlen),'C1F2').eq.1)
     &         i4head(i_bmscode)=bms_caf2
              if(index(cval(1:strlen),'Quartz').eq.1)
     &         i4head(i_bmscode)=bms_quartz
              if(index(cval(1:strlen),'SiCa').eq.1)
     &         i4head(i_bmscode)=bms_sica
              if(index(cval(1:strlen),'KBr').eq.1)
     &         i4head(i_bmscode)=bms_kbr
              if(i4head(i_bmscode).eq.0) then
                errnum=4
                write(*,'(4a)')'Error: ',path(1:lnbc(path)),
     &           ' has unexpected BMS of ',cval(1:strlen)
              elseif(verbose.ge.5) then
                write(*,*) 'BMS_cooked=',i4head(i_bmscode)
              endif
            endif
          else
            i4head(i_bmscode)=bms_caf2
            if(verbose.ge.4) then
              write(*,'(2a)') path(1:lnbc(path)),
     &         ' is missing BMS, defaulting to CaF2'
            endif
          endif

c FIXME: separate channel from detector and from AC/DC, then use these:
c      parameter (dtc_ingaas = 1)
c      parameter (dtc_si     = 2)

          call get_opus_string(luns,ipoptpar,'DTC',typ_senum,3,iendian,
     &     path,verbose,1,errnum,strlen,cval)
          if(errnum.eq.0) then
c     DG090401:  start
c     Copy DTC parameter direct to output file
            DTCstr = cval(1:strlen)
            i4head(i_dtccode)=0
c
c            if(index(cval(1:strlen),
c     &       'Si + InGaAs  DC [Internal Pos.1]').eq.1) then   !Wollongong IFS125
c              i4head(i_dtccode)=2
c            elseif(index(cval(1:strlen),
c     &       'Si + InGaAs  AC [Internal Pos.2]').eq.1) then   !Wollongong IFS125
c              i4head(i_dtccode)=1
c            elseif(index(cval(1:strlen),
c     &       'Si- Diode & Slave [Internal Pos.2]').eq.1) then
c              i4head(i_dtccode)=1
c            elseif(index(cval(1:strlen),'Si').eq.1) then
c              i4head(i_dtccode)=1
c            elseif(index(cval(1:strlen),'InGaAs').eq.1) then
c              i4head(i_dtccode)=1
c            elseif(index(cval(1:strlen),
c     &       'RT-Si Diode AC + RT-InGaAs AC [Internal Pos.2]').eq.1)then
c              i4head(i_dtccode)=1
c            elseif(index(cval(1:strlen),
c     &       'RT-Si Diode DC + RT-InGaAs DC [Internal Pos.1]').eq.1)then
c              i4head(i_dtccode)=2
c            elseif(index(cval(1:strlen),
c     &       'LN-InSb FOV=30<B0> [Internal Pos.4]').eq.1) then
c              i4head(i_dtccode)=dtc_insb
c            elseif(index(cval(1:strlen),'Bolometer').eq.1) then
c              i4head(i_dtccode)=dtc_bolom
c            else
c              errnum=4
c              write(*,'(4a)')'Error: ',path(1:lnbc(path)),
c     &         ' has unexpected DTC of ',cval(1:strlen)
c            endif
c            if((errnum.eq.0).and.(verbose.ge.5)) then
c              write(*,*) 'DTC_cooked=',i4head(i_dtccode)
c            endif
c     DG090401: end

          endif

          call test_opus_prm(luns,ipoptpar,'PGN',iendian,paramtyp)
          if(paramtyp.eq.typ_enum) then
            cval=char(0)
            call get_opus_string(luns,ipoptpar,'PGN',typ_enum,1,iendian,
     &       path,verbose,1,errnum,strlen,cval)
            if(errnum.eq.0) then
              read(cval(1:strlen),*,iostat=inpstat) i4head(i_pgn)
              if(inpstat.ne.0) then
                errnum=4
                write(*,'(3a)')'Error: ',
     &           path(1:lnbc(path)),' has incorrect PGN format'
              elseif(verbose.ge.5) then
                write(*,*) 'PGN_cooked=',i4head(i_pgn)
              endif
            endif
          else
            i4head(i_pgn)=(-1)
          endif

          call test_opus_prm(luns,ipoptpar,'HPF',iendian,paramtyp)
          if(paramtyp.eq.typ_enum) then
            cval=char(0)
            call get_opus_string(luns,ipoptpar,'HPF',typ_enum,1,iendian,
     &       path,verbose,1,errnum,strlen,cval)
            if(errnum.eq.0) then
              read(cval(1:strlen),*,iostat=inpstat) i4head(i_hpf)
              if(inpstat.ne.0) then
                errnum=4
                write(*,'(3a)')'Error: ',
     &           path(1:lnbc(path)),' has incorrect HPF format'
              elseif(verbose.ge.5) then
                write(*,*) 'HPF_cooked=',i4head(i_hpf)
              endif
            endif
          else
            i4head(i_hpf)=(-1)
          endif

          cval=char(0)
          call get_opus_string(luns,ipoptpar,'LPF',typ_enum,1,iendian,
     &     path,verbose,1,errnum,strlen,cval)
          if(errnum.eq.0) then
            read(cval(1:strlen),*,iostat=inpstat)r8head(i_lpf)
            if(inpstat.ne.0) then
              errnum=4
              write(*,'(3a)')'Error: ',
     &         path(1:lnbc(path)),' has incorrect LPF format'
            elseif(verbose.ge.5) then
              write(*,*) 'LPF_cooked=',r8head(i_lpf)
            endif
          endif

        endif        ! if((errnum.eq.0).and.(ipoptpar.ne.0))
c
c  Read from sample origin parameter block.
c
        if((errnum.eq.0).and.(iporgpar.ne.0)) then
          call get_opus_string(luns,iporgpar,'SFM',typ_string,3,iendian,
     &     path,verbose,1,errnum,strlen,cval)
          if(errnum.eq.0) then
            i4head(i_sfmcode)=0
            if(index(cval(1:strlen),'Solar').eq.1)
     &       i4head(i_sfmcode)=sfm_solar
            if(index(cval(1:strlen),'O2').gt.0)
     &       i4head(i_sfmcode)=sfm_solar
            if(index(cval(1:strlen),'Cell').gt.0)
     &       i4head(i_sfmcode)=sfm_cell
            if(index(cval(1:strlen),'cell').gt.0)
     &       i4head(i_sfmcode)=sfm_cell
            if(index(cval(1:strlen),'Tungsten lamp').eq.1)
     &       i4head(i_sfmcode)=sfm_cell
            if(index(cval(1:strlen),'Short').eq.1)
     &       i4head(i_sfmcode)=sfm_idle
            if(index(cval(1:strlen),'scripting').eq.1)
     &       i4head(i_sfmcode)=sfm_script
            if(index(cval(1:strlen),'Aerosol').eq.1)
     &       i4head(i_sfmcode)=sfm_aeros
            if(i4head(i_sfmcode).eq.0) then
              i4head(i_sfmcode)=sfm_solar
              if(verbose.ge.3) then
                write(*,'(4a)') path(1:lnbc(path)),
     &           ' has unexpected SFM of ',cval(1:strlen),
     &           ', defaulting to Solar'
              endif
            elseif(verbose.ge.5) then
              write(*,*) 'SFM_cooked=',i4head(i_sfmcode)
            endif
          endif
        endif        ! if((errnum.eq.0).and.(iporgpar.ne.0))
c
c  Read from instrument status parameter block.
c
        if((errnum.eq.0).and.(ipinstr.ne.0)) then
          call get_opus_string(luns,ipinstr,'INS',typ_string,3,iendian,
     &     path,verbose,1,errnum,strlen,cval)
          if(errnum.eq.0) then
c         DG 090602 - as for DTC above, preserve INSstring
          INSstr = cval(1:strlen)
          i4head(i_inscode)=0
c            i4head(i_inscode)=0
c            if(index(cval(1:strlen),'IFS120').eq.1)
c     &       i4head(i_inscode)=ins_120
c            if(index(cval(1:strlen),'IFS120M').eq.1)
c     &       i4head(i_inscode)=ins_120m
c            if(index(cval(1:strlen),'IFS 125').eq.1)
c     &       i4head(i_inscode)=ins_125
c            if(index(cval(1:strlen),'IFS125 HR').eq.1)
c     &       i4head(i_inscode)=ins_125
c            if(index(cval(1:strlen),'IFS 120/125 HR').eq.1)
c     &       i4head(i_inscode)=ins_125
c            if(index(cval(1:strlen),'IFS66/S').eq.1)
c     &       i4head(i_inscode)=ins_66s
c            if(index(cval(1:strlen),
c     &       'Undefined (No flange board found)').eq.1)
c     &       i4head(i_inscode)=ins_125
c            if(i4head(i_inscode).eq.0) then
c              errnum=4
c              write(*,'(4a)')'Error: ',path(1:lnbc(path)),
c     &         ' has unexpected INS of ',cval(1:strlen)
c            elseif(verbose.ge.5) then
c              write(*,*) 'INS_cooked=',i4head(i_inscode)
c            endif
c          endif
cDG110316
          elseif(errnum.eq.4)then
                INSstr='n/a'
                print *, 'INS set to "n/a"'
                errnum=0
          endif   !   if((errnum.eq.0).and.(ipinstr.ne.0))
cDG-end

          if((errnum.eq.0).and.(ipinstr.ne.0)) then
          call get_opus_i4(luns,ipinstr,'SSP',iendian,path,verbose,
     &     1,errnum,i4head(i_ssp))
cDG110316
          if(errnum.eq.4)then
             print *, 'SSP set to 1'
             i4head(i_ssp)=1
             errnum=0
          endif
          endif   !   if((errnum.eq.0).and.(ipinstr.ne.0)) then
cDG-end

cDG110316^M                                                              
        if((errnum.eq.0).and.(ipinstr.ne.0)) then
          call get_opus_i4(luns,ipinstr,'ASS',iendian,path,verbose,
     &     slicnt,errnum,ass)
            if(errnum.eq.4)then !Can't find ASS        
                ass=nss
                print *, 'Set ass=nss'
            endif
            if(ass.lt.nss)then
              print *, 'NSS =',nss, ' > ASS= ',ass
              nss=ass
            endif
         endif   ! if((errnum.eq.0).and.(ipinstr.ne.0))
cDG110316end          
          call test_opus_prm(luns,ipinstr,'SSM',iendian,paramtyp)
          if(paramtyp.eq.typ_i4) then
            call get_opus_i4(luns,ipinstr,'SSM',iendian,path,verbose,
     &       1,errnum,i4head(i_ssm))
cDG110316
              if(errnum.eq.4)then
                print *, 'SSM set to 1'
                i4head(i_ssm)=1
                errnum=0
              endif
cDG-end
          else
            i4head(i_ssm)=1
          endif

          call get_opus_r8(luns,ipinstr,'LWN',iendian,path,verbose,
     &     1,errnum,r8head(i_lwn))
cDG110316
          if(errnum.eq.4)then
                print *, 'LWN set to 15798.0138'
                r8head(i_lwn)=15798.0138
                errnum=0
            endif
cDG-end

          call test_opus_prm(luns,ipinstr,'FOC',iendian,paramtyp)
          if(paramtyp.eq.typ_r8) then
             call get_opus_r8(luns,ipinstr,'FOC',iendian,path,verbose,
     &       1,errnum,r8head(i_foc))
          elseif(errnum.ne.0) then  ! If earlier/previous error, don't rely on INS
             continue
          elseif( index(cval(1:strlen),'IFS').gt.0) then
             ii=index(cval(1:strlen),'IFS')
             if(index(cval(ii+3:strlen),'66').gt.0) then
                r8head(i_foc)=153.5d0
                if(verbose.ge.3) write(*,'(3a)') path(1:lnbc(path)),
     &          ' missing FOC. Using 215 mm based on INS of ',
     &          cval(1:strlen)
             elseif(index(cval(ii+3:strlen),'120').gt.0 .or.
     &          index(cval(ii+3:strlen),'125').gt.0 ) then
                if(index(cval(ii+6:strlen),'M').gt.0 ) then
                   r8head(i_foc)=215.0d0
                   if(verbose.ge.3) write(*,'(3a)') path(1:lnbc(path)),
     &             ' missing FOC. Using 215 mm based on INS of ',
     &             cval(1:strlen)
                else
                   r8head(i_foc)=418.0d0
                   if(verbose.ge.3) write(*,'(3a)') path(1:lnbc(path)),
     &             ' missing FOC. Using 418 mm based on INS of ',
     &             cval(1:strlen)
                endif
             else
                errnum=4
                write(*,'(3a)')'Error: unrecognised instrument ',
     &          cval(1:strlen), path(1:lnbc(path))
             endif
          else
             errnum=4
             write(*,'(3a)')'Error: unrecognised instrument ',
     &       cval(1:strlen), path(1:lnbc(path))
          endif

c          elseif((i4head(i_inscode).eq.ins_120).or.
c     &           (i4head(i_inscode).eq.ins_125)) then
c            r8head(i_foc)=418.d0
c            if(verbose.ge.3) then
c              write(*,'(3a)') path(1:lnbc(path)),' is missing FOC, ',
c     &         'using 418 mm based on INS of IFS120HR or IFS125HR'
c            endif
c          elseif(i4head(i_inscode).eq.ins_120m) then
c            r8head(i_foc)=215.d0
c            if(verbose.ge.3) then
c              write(*,'(3a)') path(1:lnbc(path)),' is missing FOC, ',
c     &         'using 215 mm based on INS of IFS120M'
c            endif
c          elseif(i4head(i_inscode).eq.ins_66s) then
c            r8head(i_foc)=153.5d0
c            if(verbose.ge.3) then
c              write(*,'(3a)') path(1:lnbc(path)),' is missing FOC, ',
c     &         'using 153.5 mm based on INS of IFS66/S'
c            endif
c          else
c            errnum=4
c            write(*,'(2a)')'Error: could not determine FOC value for ',
c     &       path(1:lnbc(path))
c          endif

          call get_opus_r8(luns,ipinstr,'DUR',iendian,path,verbose,
     &     1,errnum,r8head(i_dur))
cDG110316
          if(errnum.eq.4)then
              errnum=0
                print *, 'No information for scan time'
            endif
cDG-end


          call get_opus_i4(luns,ipinstr,'RSN',iendian,path,verbose,
     &     1,errnum,i4head(i_rsn))

          call get_opus_i4(luns,ipinstr,'PKA',iendian,path,verbose,
     &     1,errnum,i4head(i_pka))
          call get_opus_i4(luns,ipinstr,'PKL',iendian,path,verbose,
     &     1,errnum,i4head(i_pkl))

c Look for GFW and GBW, and set them to -1 if not found - DW 20130708
          call test_opus_prm(luns,ipinstr,'GFW',iendian,paramtyp)
          if(paramtyp.eq.typ_i4) then
            call get_opus_i4(luns,ipinstr,'GFW',iendian,path,verbose,
     &       1,errnum,i4head(i_gfw))
          else
            i4head(i_gfw)=-1
          endif
          call test_opus_prm(luns,ipinstr,'GBW',iendian,paramtyp)
          if(paramtyp.eq.typ_i4) then
            call get_opus_i4(luns,ipinstr,'GBW',iendian,path,verbose,
     &       1,errnum,i4head(i_gbw))
          else
            i4head(i_gbw)=-1
          endif

          call test_opus_prm(luns,ipinstr,'PRA',iendian,paramtyp)
          if(paramtyp.eq.typ_i4) then
            call get_opus_i4(luns,ipinstr,'PRA',iendian,path,verbose,
     &       1,errnum,i4head(i_pra))
          else
            i4head(i_pra)=0
          endif

          call test_opus_prm(luns,ipinstr,'PRL',iendian,paramtyp)
          if(paramtyp.eq.typ_i4) then
            call get_opus_i4(luns,ipinstr,'PRL',iendian,path,verbose,
     &       1,errnum,i4head(i_prl))
          else
            i4head(i_prl)=tpx-i4head(i_pkl)    ! FIXME ?
          endif

          call test_opus_prm(luns,ipinstr,'P2A',iendian,paramtyp)
          if(paramtyp.eq.typ_i4) then
            call get_opus_i4(luns,ipinstr,'P2A',iendian,path,verbose,
     &       1,errnum,i4head(i_p2a))
          else
c            i4head(i_p2a)=0
            i4head(i_p2a)=i4head(i_pka)
c To be able to use InGaAs-only interferograms generated by SLICE-I2S,
c we need to handle slave (channel 1) interferograms, but they contain
c PKA/PKL/PRA/PRL and not P2A/P2L/P2R/P2K - DW 20130708
          endif

          call test_opus_prm(luns,ipinstr,'P2L',iendian,paramtyp)
          if(paramtyp.eq.typ_i4) then
            call get_opus_i4(luns,ipinstr,'P2L',iendian,path,verbose,
     &       1,errnum,i4head(i_p2l))
          else
c            i4head(i_p2l)=0
            i4head(i_p2l)=i4head(i_pkl)
c DW 20130708
          endif

          call test_opus_prm(luns,ipinstr,'P2R',iendian,paramtyp)
          if(paramtyp.eq.typ_i4) then
            call get_opus_i4(luns,ipinstr,'P2R',iendian,path,verbose,
     &       1,errnum,i4head(i_p2r))
          else
c            i4head(i_p2r)=0
            i4head(i_p2r)=i4head(i_pra)
c DW 20130708
          endif

          call test_opus_prm(luns,ipinstr,'P2K',iendian,paramtyp)
          if(paramtyp.eq.typ_i4) then
            call get_opus_i4(luns,ipinstr,'P2K',iendian,path,verbose,
     &       1,errnum,i4head(i_p2k))
          else
c            i4head(i_p2k)=0
            i4head(i_p2k)=i4head(i_prl)
c DW 20130708
          endif

! The following is a pseudo-sanity check to see if the intensity
!  is sufficient to rule out a 'noise' run.
          if(i4head(i_srccode).eq.src_off) then
            if(abs(i4head(i_pka)).gt.1000.or.abs(i4head(i_pra)).gt.1000
     &       .or.abs(i4head(i_p2a)).gt.1000.or.
     &       abs(i4head(i_p2r)).gt.1000) then
              i4head(i_srccode)=src_sun
              write(*,*) 'Warning: Peak amplitude > 1000,'
              write(*,*) 'assuming Solar spectrum'
            endif
          endif

        endif        ! if((errnum.eq.0).and.(ipinstr.ne.0))

        close(unit=luns,iostat=inpstat)
        if(inpstat.ne.0) then
          errnum=4
          write(*,'(2a)')'Error: close failed on file ',path
        endif
      endif
c
c  Display date and time.  Then add stardate for the first slice or
c  time difference for all other slices.
c

ccc      if((ipoptpar.ne.0).and.(errnum.eq.0)) then
ccc  datarate=(2.d0*r8head(i_laserate)*dble(i4head(i_ssm)))/dble(i4head(i_ssp))
ccc          if(verbose.ge.4) then
ccc        write(*,*) 'Data rate=',datarate
ccc          endif
ccc      else
ccc        datarate=15000.d0             ! FIXME
ccc      endif

      if((dstatcnt.gt.0).and.(errnum.eq.0).and.(verbose.ge.3)) then
        if(slicnt.eq.1) then
          write(*,'(2a,i4.4,5(a,i2.2),a,i3.3)') 
     &     path(1:lnbc(path)),' ',
     &     iy,'-',im,'-',id,' ',hh,':',mm,':',ss,'.',ms
        elseif((timsli(slicnt)-timsli(slicnt-1)).ne.0.d0) then
          write(*,'(2a,i4.4,5(a,i2.2),a,i3.3,f8.3)') 
     &     path(1:lnbc(path)),' ',
     &     iy,'-',im,'-',id,' ',hh,':',mm,':',ss,'.',ms,
     &     timsli(slicnt)-timsli(slicnt-1)
ccc          write(*,'(2a,i4.4,5(a,i2.2),a,i3.3,f8.3,f8.3,f8.3)') 
ccc     &     path(1:lnbc(path)),' ',
ccc     &     iy,'-',im,'-',id,' ',hh,':',mm,':',ss,'.',ms,
ccc     &     timsli(slicnt)-timsli(slicnt-1),
ccc     &     dble(nptvec(slicnt))/datarate,
ccc     &     (timsli(slicnt)-timsli(slicnt-1))-
ccc     &        (dble(nptvec(slicnt))/datarate)
        endif
      endif
c
c  We're done !  Thanks for reading this far...
c
      return
      end
