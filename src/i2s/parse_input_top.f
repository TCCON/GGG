      subroutine parse_input_top(luni,mdtc,errnum,inpstat,inpath,
     & outpath,igrmode,igrmpath,phmode,phpath,dtc1,dtc2,flimit,pattern,
     & srcindic,dtcigrm,dtcspec,outfmt,delimit,minmax,stlimstd,stlimavg,
     & ylimits,xcorlim,timecorr,mscan,fftlim,sivcfreq,lsemode,
     & pco_len,pco_thresh,proclim,verbose,run_start)
c
c  Input:
c    luni           I*4    Logical Unit Number for the parameter input file
c    mdtc           I*4    Maximum number of detectors
c
c  Input/Output:
c    errnum         I*4    Error code (0=ok, <0=fatal, >0=recoverable)
c    inpstat        I*4    Value of IOSTAT return from input file read
c
c  Output:
c    inpath         C*(*)  Directory path to input igram files
c    outpath        C*(*)  Directory path for output spectrum files
c    igrmode        I*4    Interferogram saving mode (none/raw/deglitched)
c    lsemode(mdtc)  I*4    Laser sampling error type (none/slave/master/Hase/other)
c    igrmpath       C*(*)  Directory path for output interferogram files
c    phmode         C*(*)  Phase saving mode (No/yes/odd/even)
c    phpath         C*(*)  Directory path for output phase files
c    dtc1           I*4    Starting detector number to process
c    dtc2           I*4    Ending detector number to process
c    flimit         C*(*)  Name of file containing the frequency limits
c    pattern        C*(*)  Pattern for CIT file-naming convention
c    srcindic       C*(*)  List of source indicators for output files
c    dtcigrm        C*(*)  List of detector indicators for interferogram files
c    dtcspec        C*(*)  List of detector indicators for spectrum files
c    outfmt         I*4    Format selection for output file
c    delimit        C*(*)  Field delimitor for ASCII output file
c    minmax         I*4    Count of min-max pairs in ASCII output file (0=no min-max)
c    stlimavg       R*4    Limit for suntracker average intensity
c    stlimstd       R*4    Limit for suntracker intensity standard deviation
c    ylimits(2,2)   R*4    Min/Max limits on the allowed igram peaks
c    xcorlim(mdtc)  R*4    Minimum value of the ZPD cross-correlation
c    timecorr       R*4    Time correction, added to instrument time --> UT
c    mscan          I*4    Maximum number of scans analyzed in a single execution
c    fftlim(mdtc)   I*4    Maximum log-base-2 of the FFT size for each detector
c    sivcfreq(mdtc) R*8    SIV correction frequency
c    pco_thresh(mdtc) R*8  Amplitude threshold for phase correction
c    proclim        I*4    Maximum processing stage performed by the program
c    verbose        I*4    Level of verbosity for displayed messages
c    run_start      I*4    Starting run number
c
      implicit none

      integer*4 jj,ls,li,ld,
     & luni,         ! Subroutine input argument (see above)
     & mdtc,         ! Subroutine input argument (see above)
     & inpstat,      ! Subroutine input/output argument (see above)
     & errnum,       ! Subroutine input/output argument (see above)
     & igrmode,      ! Subroutine output argument (see above)
     & lsemode(mdtc), ! Subroutine output argument (see above)
     & phmode,       ! Subroutine output argument (see above)
     & dtc1,         ! Subroutine output argument (see above)
     & dtc2,         ! Subroutine output argument (see above)
     & outfmt,       ! Subroutine output argument (see above)
     & minmax,       ! Subroutine output argument (see above)
     & mscan,        ! Subroutine output argument (see above)
     & fftlim(mdtc), ! Subroutine output argument (see above)
     & pco_len(mdtc),
     & proclim,      ! Subroutine output argument (see above)
     & verbose,      ! Subroutine output argument (see above)
     & run_start,    ! Subroutine output argument (see above)
     & idtc,         ! General loop index
     & lnbc          ! Integer function Last Non-Blank Character in string

      real*4
     & timecorr,     ! Subroutine output argument (see above)
     & stlimstd,     ! suntracker relative stdev limit
     & stlimavg,     ! Suntracker avg intemsity limit
     & ylimits(2,2), ! Interferogram amplitude thresholds
     & xcorlim(mdtc) ! Subroutine output argument (see above)

      real*8
     & sivcfreq(mdtc), ! SIV Correction Frequency (cm-1)
     & pco_thresh(mdtc)  ! SIV Correction Frequency (cm-1)

      character
     & inpath*(*),   ! Subroutine output argument (see above)
     & outpath*(*),  ! Subroutine output argument (see above)
     & igrmpath*(*), ! Subroutine output argument (see above)
     & phpath*(*),   ! Subroutine output argument (see above)
     & flimit*(*),   ! Subroutine output argument (see above)
     & pattern*(*),  ! Subroutine output argument (see above)
     & srcindic*(*), ! Subroutine output argument (see above)
     & dtcigrm*(*),  ! Subroutine output argument (see above)
     & dtcspec*(*),  ! Subroutine output argument (see above)
     & delimit*(*),  ! Subroutine output argument (see above)
     & inpstr*128    ! String used to read numbers from input file

      logical*4
     & filexist      ! Keeps track of file existence
c
c  Read path to input interferograms (OPUS files or slices).
c
      call read_input_line(luni,errnum,inpstat,inpath)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)')'Error in input file couldnt read path to igrams'
         errnum=-2
      elseif(errnum.eq.0.and.inpath.eq.'0')then
         inpath=''
c      elseif(errnum.eq.0) then
c         inquire(file=inpath,exist=filexist,iostat=inpstat)
c         write(*,*) ' Here...',filexist,inpstat, inpath
c         if(inpstat.ne.0) then
c            errnum=-3
c            write(*,'(2a)')'Error: inquire failed on path ',
c     &      inpath(1:lnbc(inpath))
c         elseif(.not.filexist) then
c            errnum=-3
c            write(*,'(3a)')'Error: path to igrams ',
c     &      inpath(1:lnbc(inpath)),' does not exist'
c         endif
      endif
c
c  Read path for output files.
      call read_input_line(luni,errnum,inpstat,outpath)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)')'Error in input file: no path for spectra'
         errnum=-2
      elseif(errnum.eq.0) then
         inquire(file=outpath,exist=filexist,iostat=inpstat)
         if(inpstat.ne.0) then
            errnum=-4
            write(*,'(2a)')'Error: inquire failed on path ',
     &      outpath(1:lnbc(outpath))
c         elseif(.not.filexist) then
c            errnum=-4
c            write(*,'(3a)')'Error: path for spectra ',
c     &      outpath(1:lnbc(outpath)),' does not exist'
         endif
      endif
c
c  Read interferogram saving parameters.
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)')
     &   'Error in input file: no interferogram saving mode'
         errnum=-2
      elseif(errnum.eq.0) then
         read(inpstr,*,iostat=inpstat) igrmode
         if(inpstat.ne.0) then
            write(*,'(2a)')
     &      'Error in input file: format error in ',
     &      'interferogram saving mode'
            errnum=-2
         elseif((igrmode.lt.0).or.(igrmode.gt.2)) then
            write(*,'(2a,i0)')'Error in input file: ',
     &      'invalid interferogram saving mode of ',igrmode
            errnum=-2
         endif
      endif

      call read_input_line(luni,errnum,inpstat,igrmpath)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)')'Error in input file: no path for interferograms'
         errnum=-2
      elseif(errnum.eq.0) then
         inquire(file=igrmpath,exist=filexist,iostat=inpstat)
         if(inpstat.ne.0) then
            errnum=-5
            write(*,'(2a)')'Error: inquire failed on path ',
     &      igrmpath(1:lnbc(igrmpath))
c         elseif(.not.filexist) then
c            errnum=-5
c            write(*,'(3a)')'Error: path for interferograms ',
c     &      igrmpath(1:lnbc(igrmpath)),' does not exist'
         endif
      endif
c
c  Read phase spectrum saving parameters.
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)')
     &   'Error in input file: no phase saving mode'
         errnum=-2
      elseif(errnum.eq.0) then
         read(inpstr,*,iostat=inpstat) phmode
         if(inpstat.ne.0) then
            write(*,'(2a)')
     &      'Error in input file: format error in ',
     &      'interferogram saving mode'
            errnum=-2
         elseif((phmode.lt.0).or.(phmode.gt.1)) then
            write(*,'(2a,i0)')'Error in input file: ',
     &      'invalid interferogram saving mode of ',phmode
            errnum=-2
         endif
      endif

      call read_input_line(luni,errnum,inpstat,phpath)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)')'Error in input file: no path for phase curves'
         errnum=-2
      elseif(errnum.eq.0) then
         inquire(file=phpath,exist=filexist,iostat=inpstat)
         if(inpstat.ne.0) then
            errnum=-5
            write(*,'(2a)')'Error: inquire failed on path ',
     &      phpath(1:lnbc(phpath))
c         elseif(.not.filexist) then
c            errnum=-5
c            write(*,'(3a)')'Error: path for phase curves ',
c     &      phpath(1:lnbc(phpath)),' does not exist'
         endif
      endif
c
c  Read starting and ending detector numbers to be processed.
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)')'Error in input file: no detector specified'
         errnum=-2
      elseif(errnum.eq.0) then
         read(inpstr,*,iostat=inpstat) dtc1,dtc2
         if(inpstat.ne.0) then
            write(*,'(2a)')
     &      'Error in input file: format error in ',
     &      'detector selections'
            errnum=-2
         elseif((dtc1.lt.1).or.(dtc1.gt.mdtc)) then
            write(*,'(2(a,i0))')'Error in input file: start detector# ',
     &      dtc1,' is not between 1 and ',mdtc
            errnum=-2
         elseif((dtc2.lt.1).or.(dtc2.gt.mdtc)) then
            write(*,'(2(a,i0))')'Error in input file: end detector# ',
     &      dtc2,' is not between 1 and ',mdtc
            errnum=-2
         elseif((dtc1.gt.dtc2)) then
            write(*,'(2(a,i0))')'Error in input file: detector order: ',
     &      dtc1,' > ',dtc2
            errnum=-2
         endif
      endif
c
c  Read the full name of the frequency limit file.
      call read_input_line(luni,errnum,inpstat,flimit)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)')'Error in input file: no frequency limit file' 
         errnum=-2
      elseif(errnum.eq.0) then
         inquire(file=flimit,exist=filexist,iostat=inpstat)
         if(inpstat.ne.0) then
            errnum=-6
            write(*,'(2a)')'Error: inquire failed on file ',
     &      flimit(1:lnbc(flimit))
         elseif(.not.filexist) then
            errnum=-6
            write(*,'(3a)')'Error: frequency limit file ',
     &      flimit(1:lnbc(flimit)),' does not exist'
         endif
      endif
c
c  Read pattern for output file naming.
      call read_input_line(luni,errnum,inpstat,pattern)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)')'Error in input file: no output file pattern' 
         errnum=-2
      endif
c
c  Read source indicator.
      call read_input_line(luni,errnum,inpstat,srcindic)
      ls=lnbc(srcindic)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)')
     &   'Error in input file: no source naming pattern' 
         errnum=-2
      elseif((errnum.eq.0).and.(ls.ne.3)) then
         write(*,'(a,a,i0,a)')'Error in input file: ',
     &   'source naming pattern is of length ',ls,' instead of 3'
         errnum=-2
      endif
c
c  Read detector indicators for interferogram file naming.
      call read_input_line(luni,errnum,inpstat,dtcigrm)
      ld=lnbc(dtcigrm)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)')
     &   'Error in input file: no interferogram naming pattern' 
         errnum=-2
      elseif((errnum.eq.0).and.(ld.lt.dtc2)) then
         write(*,'(a,a,i0,a,i0)')'Error in input file.  ',
     &   'Too few detector indicators for interferograms: got ',
     &   ld,' but need ',dtc2
         errnum=-2
      endif
c
c  Read detector indicators for spectrum file naming.
      call read_input_line(luni,errnum,inpstat,dtcspec)
      ld=lnbc(dtcspec)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)')
     &   'Error in input file: no spectrum naming pattern' 
         errnum=-2
      elseif((errnum.eq.0).and.(ld.lt.dtc2)) then
         write(*,'(5a)')'Error in input file.  ',
     &   'Too few detector indicators for spectra: got ',
     &   ld,' but need ',dtc2
         errnum=-2
      endif
c
c  Read output file format selection.
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)')'Error in input file: no output format selection'
         errnum=-2
      elseif(errnum.eq.0) then
         read(inpstr,*,iostat=inpstat) outfmt
         if(inpstat.ne.0) then
            write(*,'(2a)')'Error in input file: format error in ',
     &      'output file format selection'
            errnum=-2
         elseif((outfmt.lt.0).or.(outfmt.gt.3)) then
            write(*,'(a,a,i0)')'Error in input file: ',
     &      'invalid output format of ',outfmt
            errnum=-2
         endif
      endif
c
c  Read delimitor for ASCII output files.
      call read_input_line(luni,errnum,inpstat,inpstr)
      li=lnbc(inpstr)
      if(errnum.eq.0) then
         if(inpstat.ne.0) then
            write(*,'(a)')
     &      'Error in input file: no delimitor for ASCII output files'
            errnum=-2
         elseif(li.ne.3) then
            write(*,'(2a,i6,a)')'Error in input file:  ASCII ',
     &      'delimitor selection is of length ',li,' instead of 3'
            errnum=-2
         elseif((inpstr(1:1).ne.'"').or.(inpstr(3:3).ne.'"')) then
            write(*,'(2a)')'Error in input file: ',
     &      'ASCII delimitor must be enclosed in double-quotes'
            errnum=-2
         else
            delimit=inpstr(2:2)
         endif
      endif
c
c  Read number of min-max pairs for output ASCII file.  A value of zero
c  means that output is a straight conversion of input.
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)')'Error in input file: no min-max specified'
         errnum=-2
      elseif(errnum.eq.0) then
         read(inpstr,*,iostat=inpstat) minmax
         if(inpstat.ne.0) then
            write(*,'(a)')
     &      'Error in input file: format error in min-max selection'
            errnum=-2
         elseif(minmax.lt.0) then
            write(*,'(2a,i0)')'Error in input file: ',
     &      'invalid minmax specification of ',minmax
            errnum=-2
         endif
      endif
c
c  Read suntracker total intensity limits (STD and AVG).
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)')
     &   'Error in input file: no suntracker intensity limits'
         errnum=-2
      elseif(errnum.eq.0) then
         read(inpstr,*,iostat=inpstat) stlimstd,stlimavg
         if(inpstat.ne.0) then
            write(*,'(2a)')
     &      'Error in input file: format error in ',
     &      'suntracker intensity limits'
            errnum=-2
         elseif((stlimstd.lt.0.0).or.(stlimavg.lt.0.0)) then
            write(*,'(2a)')'Error in input file: ',
     &      'suntracker intensity limits must be positive'
            errnum=-2
         endif
      endif

c
c  Read min and max allowed values of the interferogram.
      do jj=1,2  ! Min and Max
         call read_input_line(luni,errnum,inpstat,inpstr)
         if(inpstat.ne.0) then
            write(*,'(a)')
     &      'Error in input file: no igram min/max values '
            errnum=-2
         else
            read(inpstr,*,iostat=inpstat)(ylimits(idtc,jj),idtc=1,dtc2)
            if(inpstat.ne.0) then
               write(*,'(2a)')
     &         'Error in input file: format error in ',
     &         'min/max values for the igram'
               errnum=-2
            endif
         endif
      end do   !  do jj=1,2  ! Min and Max
c
c  Read minimum value of the cross-correlation between runs.
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)')
     &   'Error in input file: no values for minimum cross-correlation'
         errnum=-2
      elseif(errnum.eq.0) then
         read(inpstr,*,iostat=inpstat)(xcorlim(idtc),idtc=1,dtc2)
         if(inpstat.ne.0) then
            write(*,'(2a)')
     &      'Error in input file: format error in ',
     &      'values for minimum cross-correlation'
            errnum=-2
         else
            do idtc=1,dtc2
               if((errnum.eq.0).and.(xcorlim(idtc).lt.0.0)) then
                  write(*,'(2a)')'Error in input file: ',
     &            'cross-correlation limits must be positive'
                  errnum=-2
               endif
            enddo
         endif
      endif

c
c  Read the time correction, in floating-point hours.
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)') 'Error in input file: no time correction'
         errnum=-2
      elseif(errnum.eq.0) then
         read(inpstr,*,iostat=inpstat) timecorr
         if(inpstat.ne.0) then
            write(*,'(a)')
     &      'Error in input file: format error in time correction'
            errnum=-2
         else
            if((errnum.eq.0).and.
     &      ((timecorr.lt.(-24.0)).or.(timecorr.gt.24.0))) then
               write(*,'(2a)')'Error in input file: ',
     &         'time correction must be between -24 and +24'
               errnum=-2
            endif
         endif
      endif
c
c  Read maximum number of scans to analyze.
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)')'Error in input file: no max number of scans'
         errnum=-2
      elseif(errnum.eq.0) then
         read(inpstr,*,iostat=inpstat) mscan
         if(inpstat.ne.0) then
            write(*,'(2a)')'Error in input file: format error in ',
     &      'maximum number of scans'
            errnum=-2
         elseif(mscan.lt.0) then
            write(*,'(2a,i0)')'Error in input file: ',
     &      'invalid maximum number of scans of ',mscan
            errnum=-2
         endif
      endif
c
c  Read maximum long-side igram size.
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)')
     &   'Error in input file: no values for maximum FFT size'
         errnum=-2
      elseif(errnum.eq.0) then
         read(inpstr,*,iostat=inpstat)(fftlim(idtc),idtc=1,dtc2) 
         if(inpstat.ne.0) then
            write(*,'(a)')'Error in input file: format error: Mlong'
            errnum=-2
         else
            do idtc=1,dtc2
               if((errnum.eq.0).and.(fftlim(idtc).lt.0)) then
                  write(*,'(a)')'Error in input file: Mlong must be +ve'
                  errnum=-2
               endif
            enddo
         endif
      endif
c
c  Read SIV_correction frequencies
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)')
     &   'Error in input file: no values for SIVCF'
         errnum=-2
      elseif(errnum.eq.0) then
         read(inpstr,*,iostat=inpstat)(sivcfreq(idtc),idtc=1,dtc2) 
         if(inpstat.ne.0) then
            write(*,'(2a)')'Error in input file: format error in ',
     &      'values of SIV correction frequencies'
            errnum=-2
         else
            do idtc=1,dtc2
               if((errnum.eq.0).and.(sivcfreq(idtc).lt.0)) then
                  write(*,'(2a)')'Error in input file: ',
     &           'values of SIV correction frequencies must be positive'
                  errnum=-2
               endif
            enddo
         endif
      endif
c
c  Read LSE type
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)')
     &   'Error in input file: no laser sampling error type'
         errnum=-2
      elseif(errnum.eq.0) then
         read(inpstr,*,iostat=inpstat)(lsemode(idtc),idtc=1,dtc2)
         if(inpstat.ne.0) then
            write(*,'(2a)')
     &      'Error in input file: format error in ',
     &      'laser sampling error type'
            errnum=-2
         else
            do idtc=1,dtc2
               if((errnum.eq.0).and.
     &         ((lsemode(idtc).lt.0).or.(lsemode(idtc).gt.4))) then
                  write(*,'(2a,i0)')'Error in input file: ',
     &           'invalid laser sampling error type of ',lsemode(idtc)
                  errnum=-2
               endif
            enddo
         endif
      endif

c
c  Read lengths of the PCO
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)') 'Error in input file: no PCO thresholds'
         errnum=-2
      elseif(errnum.eq.0) then
         read(inpstr,*,iostat=inpstat)(pco_len(idtc),idtc=1,dtc2) 
         if(inpstat.ne.0) then
            write(*,'(2a)')
     &     'Error in input file: format error in PCO lengths'
            errnum=-2
         else
            do idtc=1,dtc2
               if((errnum.eq.0).and.(pco_len(idtc).lt.0)) then
                  write(*,'(2a)')'Error in input file: ',
     &           'values of PCO lengths must be positive'
                  errnum=-2
               endif
            enddo
         endif
      endif
c
c  Read spectral magnitude thresholds for the PCO
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)') 'Error in input file: no PCO thresholds'
         errnum=-2
      elseif(errnum.eq.0) then
         read(inpstr,*,iostat=inpstat)(pco_thresh(idtc),idtc=1,dtc2) 
         if(inpstat.ne.0) then
            write(*,'(2a)')
     &      'Error in input file: format error in PCO thresholds'
            errnum=-2
         else
            do idtc=1,dtc2
               if((errnum.eq.0).and.(pco_thresh(idtc).lt.0)) then
                  write(*,'(2a)')'Error in input file: ',
     &            'values of PCO thresholds must be positive'
                  errnum=-2
               endif
            enddo
         endif
      endif
c
c  Read maximum level of processing.
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)')
     &   'Error in input file: no maximum level of processing'
         errnum=-2
      elseif(errnum.eq.0) then
         read(inpstr,*,iostat=inpstat) proclim
         if(inpstat.ne.0) then
            write(*,'(2a)')
     &     'Error in input file: format error in ',
     &     'maximum level of processing'
            errnum=-2
         elseif((proclim.lt.0).or.(proclim.gt.4)) then
            write(*,'(2a,i0)')'Error in input file: ',
     &      'invalid maximum level of processing of ',proclim
            errnum=-2
         endif
      endif
c
c  Read verbosity level.
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)')'Error in input file: no verbosity specified'
         errnum=-2
      elseif(errnum.eq.0) then
         read(inpstr,*,iostat=inpstat) verbose
         if(inpstat.ne.0) then
            write(*,'(2a)')
     &      'Error in input file: format error in ',
     &      'verbosity specification'
            errnum=-2
         elseif((verbose.lt.1).or.(verbose.gt.5)) then
            write(*,'(2a,i0)')'Error in input file: ',
     &      'invalid verbosity specification of ',verbose
            errnum=-2
         endif
      endif
c
c  Read verbosity level.
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)')'Error in input file: no run_start specified'
         errnum=-2
      elseif(errnum.eq.0) then
         read(inpstr,*,iostat=inpstat) run_start
         if(inpstat.ne.0) then
            write(*,'(2a)')
     &      'Error in input file: format error in ',
     &      'verbosity specification'
            errnum=-2
         elseif(run_start.lt.1) then
            write(*,'(2a,i0)')'Error in input file: ',
     &      'invalid run_start specification of ',verbose
            errnum=-2
         endif
      endif
c
      return
      end
