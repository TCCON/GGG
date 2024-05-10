      subroutine get_igram_run_parameters(inpath,fpfn,
     & catslice,runno,verbose,mns,msl,mip,mch,mi4,mr8,
     & chan1,chan2,errnum,nptvec,bpdata,nss,tpx,timvec,
     & timsli,runsta,runend,i4head,r8head,DTCstr,INSstr)
c
c  Input:
c    inpath        C*(*)  Directory path to input slice files
c    fpfn          C*(*)  Name of input OPUS file
c    catslice      I*4    Slice number from catalog, used in slice file name
c    runno         I*4    Run number increasing throughout the day
c    verbose       I*4    Level of verbosity for displayed messages
c    mns           I*4    Maximum number of interferograms per scan set (max NSS)
c    msl           I*4    Maximum number of interferogram slices per scan set
c    mip           I*4    Maximum number of input points
c    mch           I*4    Maximum number of data channels
c    mi4           I*4    Maximum number of I*4 items in file header
c    mr8           I*4    Maximum number of R*8 items in file header
c    chan1         I*4    Starting channel number to process
c    chan2         I*4    Ending channel number to process
c
c  Input/output:
c    errnum        I*4    Error code (0=ok, <0=fatal, >0=recoverable)
c
c  Output:
c    nptvec(msl)     I*4  Number of PoinTs, from OPUS header
c    bpdata(msl,mch) I*4  Byte pointers into the data blocks of slices
c    nss             I*4  Number of Sample Scans in this set of slices
c    tpx             I*4  Number of points in one FWD or REV scan
c    timvec(mns)     R*8  Time vector contains one entry for each scan
c    timsli(msl)     R*8  Time vector contains one entry for each slice
c    runsta(mns)     I*4  Run starting slice number ( > 0) or run error ( < 0)
c    runend(mns)     I*4  Run ending slice number
c    i4head(mi4)     I*4  Vector to hold the I*4 header items
c    r8head(mr8)     R*8  Vector to hold the R*8 header items
c    DTCstr          C*40 Variable to hold detector description
c    INSstr          C*40 Variable to hold detector description
c
      implicit none

      integer*4 idum,
     & catslice,        ! Subroutine input argument (see above)
     & runno,        ! Subroutine input argument (see above)
     & verbose,      ! Subroutine input argument (see above)
     & mns,          ! Subroutine input argument (see above)
     & msl,          ! Subroutine input argument (see above)
     & mip,          ! Subroutine input argument (see above)
     & mch,          ! Subroutine input argument (see above)
     & mi4,          ! Subroutine input argument (see above)
     & mr8,          ! Subroutine input argument (see above)
     & chan1,        ! Subroutine input argument (see above)
     & chan2,        ! Subroutine input argument (see above)
     & errnum,       ! Subroutine input/output argument (see above)
     & nptvec(msl),  ! Subroutine output argument (see above)
     & bpdata(msl,mch),! Subroutine output argument (see above)
     & nss,          ! Subroutine output argument (see above)
     & tpx,          ! Subroutine output argument (see above)
     & runsta(mns),  ! Subroutine output argument (see above)
     & runend(mns),  ! Subroutine output argument (see above)
     & i4head(mi4),  ! Subroutine output argument (see above)
     & slicnt,       ! Slice count, number of elements of npt and bpdata vectors
     & lnbc,         ! Integer function Last Non-Blank Character in string
c     & fnbc,        ! Integer function First Non-Blank Character in string
     & lf,
     & slicind,      ! Slice index
     & kss,          ! Count of scans
     & counter,      ! Count of data points
     & kk,
     & errmem,       ! Value of errnum at start of scan splitting
     & ijk,       ! General loop index
     & jd,j1         ! Julian Days

      real*8
     & timvec(mns),  ! Subroutine output argument (see above)
     & timsli(msl),  ! Subroutine output argument (see above)
     & r8head(mr8),  ! Subroutine output argument (see above)
     & firminst      ! Time of the install of the "1.300 Apr 14 2004" firmware

      character
     & inpath*(*),   ! Subroutine input argument (see above)
     & fpfn*(*),     ! Subroutine input argument (see above)
     & filename*256, ! Full file name for data slices
     & DTCstr*(*),   ! Character variable to hold DTC description^M
     & INSstr*(*)    ! Character variable to hold INS description^M

      logical*4
     & filexist,     ! Keeps track of file existence
     & inrun         ! True while inside NSS interferometer runs

      include 'header_indices.inc'

      idum=bms_caf2  ! Avoid compiler warnings (unused)
      idum=bms_kbr   ! Avoid compiler warnings (unused)
      idum=bms_quartz! Avoid compiler warnings (unused)
      idum=bms_sica  ! Avoid compiler warnings (unused)
      idum=i_aptval  ! Avoid compiler warnings (unused)
      idum=i_bbw     ! Avoid compiler warnings (unused)
      idum=i_bfw     ! Avoid compiler warnings (unused)
      idum=i_aqmcode ! Avoid compiler warnings (unused)
      idum=i_bmscode ! Avoid compiler warnings (unused)
      idum=i_dtccode ! Avoid compiler warnings (unused)
      idum=i_inscode ! Avoid compiler warnings (unused)
      idum=i_laserate! Avoid compiler warnings (unused)
      idum=i_dur     ! Avoid compiler warnings (unused)
      idum=i_foc     ! Avoid compiler warnings (unused)
      idum=i_hpf     ! Avoid compiler warnings (unused)
      idum=i_gfw     ! Avoid compiler warnings (unused)
      idum=i_gbw     ! Avoid compiler warnings (unused)
      idum=i_lpf     ! Avoid compiler warnings (unused)
      idum=i_mvd     ! Avoid compiler warnings (unused)
      idum=i_p2a     ! Avoid compiler warnings (unused)
      idum=i_p2k     ! Avoid compiler warnings (unused)
      idum=i_p2l     ! Avoid compiler warnings (unused)
      idum=i_p2r     ! Avoid compiler warnings (unused)
      idum=i_pgn     ! Avoid compiler warnings (unused)
      idum=i_pka     ! Avoid compiler warnings (unused)
      idum=i_pra     ! Avoid compiler warnings (unused)
      idum=i_prl     ! Avoid compiler warnings (unused)
      idum=i_pkl     ! Avoid compiler warnings (unused)
      idum=i_rsn     ! Avoid compiler warnings (unused)
      idum=i_sgna    ! Avoid compiler warnings (unused)
      idum=i_sgnb    ! Avoid compiler warnings (unused)
      idum=i_srccode ! Avoid compiler warnings (unused)
      idum=i_ssm     ! Avoid compiler warnings (unused)
      idum=i_ssp     ! Avoid compiler warnings (unused)
      idum=i_zff     ! Avoid compiler warnings (unused)
      idum=src_sun   ! Avoid compiler warnings (unused)
      idum=src_off   ! Avoid compiler warnings (unused)
      idum=src_nir   ! Avoid compiler warnings (unused)
      idum=src_mir   ! Avoid compiler warnings (unused)
      idum=sfm_solar ! Avoid compiler warnings (unused)
      idum=sfm_script! Avoid compiler warnings (unused)
      idum=sfm_cell  ! Avoid compiler warnings (unused)
      idum=sfm_aeros ! Avoid compiler warnings (unused)
      idum=aqm_sd    ! Avoid compiler warnings (unused)
      idum=aqm_sf    ! Avoid compiler warnings (unused)
      idum=aqm_sn    ! Avoid compiler warnings (unused)

c  Write a blank line to the screen to separate scan sets.
      if(verbose.ge.3)  write(*,'("")')

c  Initialize variables.
      inrun=.true.
      slicnt=0
      i4head(i_sfmcode)=0  ! Must be initialized: it is used to reject scans.
 
c  Note: slicnt is incremented inside get_opusigram_params
      do while( inrun .and. errnum.eq.0 )
 
c  Check that the next OPUS slice/file exists.
         write(filename,'(2a)')inpath(:lnbc(inpath)),fpfn(:lnbc(fpfn))
         lf=lnbc(filename)
         if(catslice.gt.0)
     &   write(filename(lf+1:),'(i0,a2)')catslice+slicnt,'.0'
         inquire(file=filename,exist=filexist)
         if(.not.filexist) then
            errnum=slicnt+1
            if(slicnt.eq.0)write(*,'(2a,i0,a)')'Error:OPUS file/slice ',
     &      filename(1:lf),catslice+slicnt,' is missing or incomplete'
            exit
         endif

c  Read full set of parameters from this OPUS file/slice.
c  Note: 'get_opusigram_params' sets 'inrun' to false at end of scan set.
c         write(*,*)'Calling get_opusigram_param',slicnt
         call get_opusigram_params(filename,verbose,msl,mch,mi4,mr8,
     &   chan1,chan2,errnum,slicnt,inrun,nptvec,bpdata,timsli,nss,tpx,
     &   i4head,r8head,DTCstr,INSstr)
c         write(*,*)'Called get_opusigram_param, line 168:',errnum,slicnt

      enddo       ! while( inrun .and. errnum.eq.0 )
 
c  If OPUS files/slices cover more than one scan, divide TPX by two
c  to separate FWD from REV runs.
      if( errnum.eq.0 .and. nss.gt.1 .and. slicnt.gt.1 ) then
         if(mod(tpx,2).ne.0) then
            errnum=2
            write(*,'(3a,i0)')'Error: ',filename(1:lnbc(filename)),
     &      ' has odd TPX of ',tpx
         else
            tpx=tpx/2
         endif
      endif
c      write(*,*) 'line 183  tpx=',tpx
 
c  Make sure TPX does not exceed MIP, the maximum number of input points.
      if( errnum.eq.0 .and. tpx.gt.mip ) then
         errnum=2
         write(*,'(a,i0)')'Error: increase parameter MIP to ',tpx
      endif
 
c  Fix laser wavenumber: IFS125-1 instrument returns an incorrect value
      if( errnum.eq.0 .and.
     &    dabs(r8head(i_lwn)-15798.1d0) .lt. 0.01d0 ) then
         r8head(i_lwn)=15798.0138d0
         if(verbose.ge.3) write(*,'(a,f11.4)')
     &   'Replacing laser wavenumber of 15798.1 with',r8head(i_lwn)
      endif
 
c  Initialize scan-splitting variables.
      counter=0
      kss=0
      kk=0
      errmem=errnum
      call julian(2004,8,11,jd)    ! Firmware "upgraded" on 11-Aug-2004
      call julian(2000,1,1,j1)
      firminst=dble((((jd-j1)*24.0d0)*60)*60)  ! this could be a parameter
 
c  Loop over the slices in this scan set to separate the individual scans.
      do slicind=1,slicnt
         if(errnum.eq.0) then
c         do while (errnum.eq.0 .and. slicind.le.slicnt)
            counter=counter+nptvec(slicind)

c  Added following 4 lines of code as per JFB email of 2018-07-11
            if((counter.eq.(2*tpx)).and.(slicnt.eq.1)) then
              write(*,*)'SINGLE SLICE MULTI SCAN RUN - NO WARRANTEE'
              counter=tpx
            endif

            if( counter.eq.tpx .and. kss.eq.mns ) then
               errnum=5
               write(*,'(a)') 'Error: increase parameter MNS'
            elseif(counter.eq.tpx) then
               if(verbose.ge.3) then
                  write(*,'(a4,i0,a26,i0)') 'Run ',runno+kss,
     &            ' found, starting at slice ',catslice+kk
                  write(*,*)'+++++++++++++++++++++'
               endif
               kss=kss+1
               runsta(kss)=0   ! Indicates no error found so far
 
c  Reject idle scans.
               if( runsta(kss).eq.0 .and.
     &         i4head(i_sfmcode).eq.sfm_idle ) then
                  runsta(kss)=-100
                  write(*,'(3a)') 'Reject: this run is an idle scan'
               endif
 
c  Reject scans affected by download stress: these have incorrect time.
               do ijk=2,kk+2
                  if(runsta(kss).eq.0) then
                     if(((timsli(ijk+1)-timsli(ijk)).gt.23.d0)
     &               .and.((timsli(ijk)-timsli(ijk-1)).ne.0.d0)) then
                        runsta(kss)=-103
                        write(*,'(3a)') 'Reject: ',filename,'
     &                  affected by download stress '
                     endif
                  endif
               enddo  ! do ijk=2,kk+2
 
c  Set starting+ending slices and time.
c  Apr 23, 2001, fix for scans with only 1 slice containing DAT/TIM info
c               write(*,*) 'slicind = ',slicind
               if(slicind.eq.1) then
                  timvec(kss)=timsli(slicind)
                  runsta(kss)=catslice+kk         ! DW 2018-07-09
                  runend(kss)=catslice+slicind-1  ! DW 2018-07-09
               elseif(runsta(kss).eq.0) then
                  runsta(kss)=catslice+kk
                  runend(kss)=catslice+slicind-1
c                 write(*,*)'kss,runsta(kss),runend(kss) = ',kss,
c    &            runsta(kss),runend(kss)
                  if( runend(kss)-runsta(kss) .lt. 2) then
                     timvec(kss)=timsli(kk+1)  ! Low resolution scan
                  elseif(timsli(slicind).lt.firminst) then ! B4 firmware upgrade
                     timvec(kss)=(2.d0*timsli(kk+2))-timsli(kk+3)
                  else                                     ! After firmware upgrade
                     timvec(kss)=(3.d0*timsli(kk+2))-(2.d0*timsli(kk+3))
                  endif
               endif

               counter=0
               kk=slicind

c  If the point count exceeds the expected value of TPX, the run is
c  not ending on a slice boundary.  This is an error.
            elseif(counter.gt.tpx) then
               errnum=5
               write(*,'(a14,i0,a8,i0,a38,i0,a17,i0)')'Error in run ',
     &         runno+kss,': slice ',catslice+slicind-1,
     &         ' goes past the end of scan.  Expected ',
     &         tpx,' points, but got ',counter
            endif  !  if( counter.eq.tpx .and. kss.eq.mns ) then
         endif
      enddo  ! do while ( errnum.eq.0 .and. slicind.le.slicnt )

c  At end-of-set check that all data was processed.  If not, the
c  last scan ended prematurely.
      if( errnum.eq.0 .and. counter.ne.0 ) then
         errnum=5
         write(*,'(a13,i0,a8,i0,a42,i0,a17,i0)')'Error in run ',
     &   runno+kss,': catslice ',catslice+slicnt,
     &   ' caused early scan termination.  Expected ',
     &   tpx,' points, but got ',counter
      endif
 
c  If there was an error, flag undetected runs as invalid.
      do while (kss.lt.nss)
         kss=kss+1
c         timvec(kss)=0.d0 !  Causes name problem with opus-i2s-example
         if(kss.gt.1) timvec(kss)=timvec(kss-1)  ! GCT kludge
         runsta(kss)=(-errnum)
      enddo
 
c  Restore error state to its value before scan splitting so that runs
c  before the error condition (if any) can be processed.
      errnum=errmem

      return
      end
