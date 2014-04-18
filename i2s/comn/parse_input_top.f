      subroutine parse_input_top(luni,mdtc,errnum,inpstat,inpath,
     & outpath,igrmode,igrmpath,phmode,phpath,dtc1,dtc2,flimit,pattern,
     & srcindic,dtcigrm,dtcspec,outfmt,delimit,minmax,stlimstd,stlimavg,
     & zpdlim,xcorlim,timecorr,scanlim,fftlim,sivcfreq,lsemode,
     & pco_len,pco_thresh,proclim,verbose)
c
c  Input:
c    luni          I*4    Logical Unit Number for the parameter input file
c    mdtc          I*4    Maximum number of detectors
c
c  Input/Output:
c    errnum        I*4    Error code (0=ok, <0=fatal, >0=recoverable)
c    inpstat       I*4    Value of IOSTAT return from input file read
c
c  Output:
c    inpath        C*(*)  Directory path to input igram files
c    outpath       C*(*)  Directory path for output spectrum files
c    igrmode       I*4    Interferogram saving mode (none/raw/deglitched)
c    lsemode(mdtc) I*4    Laser sampling error type (none/slave/master/Hase/other)
c    igrmpath      C*(*)  Directory path for output interferogram files
c    phmode        C*(*)  Phase saving mode (No/yes/odd/even)
c    phpath        C*(*)  Directory path for output phase files
c    dtc1          I*4    Starting detector number to process
c    dtc2          I*4    Ending detector number to process
c    flimit        C*(*)  Name of file containing the frequency limits
c    pattern       C*(*)  Pattern for CIT file-naming convention
c    srcindic      C*(*)  List of source indicators for output files
c    dtcigrm       C*(*)  List of detector indicators for interferogram files
c    dtcspec       C*(*)  List of detector indicators for spectrum files
c    outfmt        I*4    Format selection for output file
c    delimit       C*(*)  Field delimitor for ASCII output file
c    minmax        I*4    Count of min-max pairs in ASCII output file (0=no min-max)
c    stlimavg      R*4    Limit for suntracker average intensity
c    stlimstd      R*4    Limit for suntracker intensity standard deviation
c    zpdlim(mdtc)  R*4    Minimum value of the ZPD peaks
c    xcorlim(mdtc) R*4    Minimum value of the ZPD cross-correlation
c    timecorr      R*4    Time correction, added to instrument time --> UT
c    scanlim       I*4    Maximum number of scans analyzed in a single execution
c    fftlim(mdtc)  I*4    Maximum log-base-2 of the FFT size for each detector
c    sivcfreq(mdtc)  R*8    SIV correction frequency
c    pco_thresh(mdtc) R*8 
c    proclim       I*4    Maximum processing stage performed by the program
c    verbose       I*4    Level of verbosity for displayed messages
c
      implicit none

      integer*4
     & luni,        ! Subroutine input argument (see above)
     & mdtc,        ! Subroutine input argument (see above)
     & inpstat,     ! Subroutine input/output argument (see above)
     & errnum,      ! Subroutine input/output argument (see above)
     & igrmode,     ! Subroutine output argument (see above)
     & lsemode(mdtc), ! Subroutine output argument (see above)
     & phmode,      ! Subroutine output argument (see above)
     & dtc1,        ! Subroutine output argument (see above)
     & dtc2,        ! Subroutine output argument (see above)
     & outfmt,      ! Subroutine output argument (see above)
     & minmax,      ! Subroutine output argument (see above)
     & scanlim,     ! Subroutine output argument (see above)
     & fftlim(mdtc),! Subroutine output argument (see above)
     & pco_len(mdtc),
     & proclim,     ! Subroutine output argument (see above)
     & verbose,     ! Subroutine output argument (see above)
     & indexa,      ! General loop index
     & lnbc,        ! Integer function Last Non-Blank Character in string
     & fnbc         ! Integer function First Non-Blank Character in string

      real*4
     & timecorr,    ! Subroutine output argument (see above)
     & stlimstd,    ! suntracker relative stdev limit
     & stlimavg,    ! Suntracker avg intemsity limit
     & zpdlim(mdtc),! Interferogram amplitude threshold
     & xcorlim(mdtc)! Subroutine output argument (see above)

      real*8
     & sivcfreq(mdtc), ! SIV Correction Frequency (cm-1)
     & pco_thresh(mdtc)  ! SIV Correction Frequency (cm-1)

      character
     & inpath*(*),  ! Subroutine output argument (see above)
     & outpath*(*), ! Subroutine output argument (see above)
     & igrmpath*(*),! Subroutine output argument (see above)
     & phpath*(*),  ! Subroutine output argument (see above)
     & flimit*(*),  ! Subroutine output argument (see above)
     & pattern*(*), ! Subroutine output argument (see above)
     & srcindic*(*),! Subroutine output argument (see above)
     & dtcigrm*(*), ! Subroutine output argument (see above)
     & dtcspec*(*), ! Subroutine output argument (see above)
     & delimit*(*), ! Subroutine output argument (see above)
     & inpstr*99,   ! String used to read numbers from input file
     & stringa*11,  ! String used to format integer display
     & stringb*11   ! String used to format integer display

      logical*4
     & filexist     ! Keeps track of file existence
c
c  Read path to input interferograms (OPUS files or slices).
c
      call read_input_line(luni,errnum,inpstat,inpath)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
      write(*,'(a)')'Error in input file: could not read path to igrams'
        errnum=-2
      elseif(errnum.eq.0.and.inpath.eq.'0')then
        inpath=''
c      elseif(errnum.eq.0) then
c        inquire(file=inpath,exist=filexist,iostat=inpstat)
c        write(*,*) ' Here...',filexist,inpstat, inpath
c        if(inpstat.ne.0) then
c          errnum=-3
c          write(*,'(2a)')'Error: inquire failed on path ',
c     &     inpath(1:lnbc(inpath))
c        elseif(.not.filexist) then
c          errnum=-3
c          write(*,'(3a)')'Error: path to igrams ',
c     &     inpath(1:lnbc(inpath)),' does not exist'
c        endif
      endif
c
c  Read path for output files.
c
      call read_input_line(luni,errnum,inpstat,outpath)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
        write(*,'(a)')'Error in input file: no path for spectra'
        errnum=-2
      elseif(errnum.eq.0) then
        inquire(file=outpath,exist=filexist,iostat=inpstat)
        if(inpstat.ne.0) then
          errnum=-4
          write(*,'(2a)')'Error: inquire failed on path ',
     &     outpath(1:lnbc(outpath))
c        elseif(.not.filexist) then
c          errnum=-4
c          write(*,'(3a)')'Error: path for spectra ',
c     &     outpath(1:lnbc(outpath)),' does not exist'
        endif
      endif
c
c  Read interferogram saving parameters.
c
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
            write(stringa,'(i11)') igrmode
            write(*,'(3a)')'Error in input file: ',
     &      'invalid interferogram saving mode of ',
     &      stringa(fnbc(stringa):11)
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
c
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
            write(stringa,'(i11)') phmode
            write(*,'(3a)')'Error in input file: ',
     &      'invalid interferogram saving mode of ',
     &      stringa(fnbc(stringa):11)
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
c
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
        write(*,'(a)')'Error in input file: no detector specified'
        errnum=-2
      elseif(errnum.eq.0) then
        read(inpstr,*,iostat=inpstat) dtc1,dtc2
        if(inpstat.ne.0) then
          write(*,'(2a)')
     &     'Error in input file: format error in ',
     &     'detector selections'
          errnum=-2
        elseif((dtc1.lt.1).or.(dtc1.gt.mdtc)) then
          write(stringa,'(i11)') dtc1
          write(stringb,'(i11)') mdtc
          write(*,'(4a)')'Error in input file: starting detector of ',
     &     stringa(fnbc(stringa):11),
     &     ' is not between 1 and ',stringb(fnbc(stringb):11)
          errnum=-2
        elseif((dtc2.lt.1).or.(dtc2.gt.mdtc)) then
          write(stringa,'(i11)') dtc2
          write(stringb,'(i11)') mdtc
          write(*,'(4a)')'Error in input file: ending detector of ',
     &     stringa(fnbc(stringa):11),
     &     ' is not between 1 and ',stringb(fnbc(stringb):11)
          errnum=-2
        elseif((dtc1.gt.dtc2)) then
          write(stringa,'(i11)') dtc1
          write(stringb,'(i11)') dtc2
          write(*,'(5a)')'Error in input file.  ',
     &     'Incorrect detector order: ',
     &     stringa(fnbc(stringa):11),
     &     ' is after ',stringb(fnbc(stringb):11)
          errnum=-2
        endif
      endif
c
c  Read the full name of the frequency limit file.
c
      call read_input_line(luni,errnum,inpstat,flimit)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
        write(*,'(a)')'Error in input file: no frequency limit file' 
        errnum=-2
      elseif(errnum.eq.0) then
        inquire(file=flimit,exist=filexist,iostat=inpstat)
        if(inpstat.ne.0) then
          errnum=-6
          write(*,'(2a)')'Error: inquire failed on file ',
     &     flimit(1:lnbc(flimit))
        elseif(.not.filexist) then
          errnum=-6
          write(*,'(3a)')'Error: frequency limit file ',
     &     flimit(1:lnbc(flimit)),' does not exist'
        endif
      endif
c
c  Read pattern for output file naming.
c
      call read_input_line(luni,errnum,inpstat,pattern)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
        write(*,'(a)')'Error in input file: no output file pattern' 
        errnum=-2
      endif
c
c  Read source indicator.
c
      call read_input_line(luni,errnum,inpstat,srcindic)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
        write(*,'(a)')
     &   'Error in input file: no source naming pattern' 
        errnum=-2
      elseif((errnum.eq.0).and.(lnbc(srcindic).ne.3)) then
        write(stringa,'(i11)') lnbc(srcindic)
        write(*,'(4a)')'Error in input file: ',
     &   'source naming pattern is of length ',
     &   stringa(fnbc(stringa):11),' instead of 3'
        errnum=-2
      endif
c
c  Read detector indicators for interferogram file naming.
c
      call read_input_line(luni,errnum,inpstat,dtcigrm)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
        write(*,'(a)')
     &   'Error in input file: no interferogram naming pattern' 
        errnum=-2
      elseif((errnum.eq.0).and.(lnbc(dtcigrm).lt.dtc2)) then
        write(stringa,'(i11)') lnbc(dtcigrm)
        write(stringb,'(i11)') dtc2
        write(*,'(5a)')'Error in input file.  ',
     &   'Too few detector indicators for interferograms: got ',
     &   stringa(fnbc(stringa):11),
     &   ' but need ',stringb(fnbc(stringb):11)
        errnum=-2
      endif
c
c  Read detector indicators for spectrum file naming.
c
      call read_input_line(luni,errnum,inpstat,dtcspec)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
        write(*,'(a)')
     &   'Error in input file: no spectrum naming pattern' 
        errnum=-2
      elseif((errnum.eq.0).and.(lnbc(dtcspec).lt.dtc2)) then
        write(stringa,'(i11)') lnbc(dtcspec)
        write(stringb,'(i11)') dtc2
        write(*,'(5a)')'Error in input file.  ',
     &   'Too few detector indicators for spectra: got ',
     &   stringa(fnbc(stringa):11),
     &   ' but need ',stringb(fnbc(stringb):11)
        errnum=-2
      endif
c
c  Read output file format selection.
c
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
        write(*,'(a)')'Error in input file: no output format selection'
        errnum=-2
      elseif(errnum.eq.0) then
        read(inpstr,*,iostat=inpstat) outfmt
        if(inpstat.ne.0) then
          write(*,'(2a)')
     &     'Error in input file: format error in ',
     &     'output file format selection'
          errnum=-2
        elseif((outfmt.lt.0).or.(outfmt.gt.3)) then
          write(stringa,'(i11)') outfmt
          write(*,'(3a)')'Error in input file: ',
     &     'invalid output format of ',
     &     stringa(fnbc(stringa):11)
          errnum=-2
        endif
      endif
c
c  Read delimitor for ASCII output files.
c
      call read_input_line(luni,errnum,inpstat,inpstr)
      if(errnum.eq.0) then
        if(inpstat.ne.0) then
          write(*,'(a)')
     &     'Error in input file: no delimitor for ASCII output files'
          errnum=-2
        elseif(lnbc(inpstr).ne.3) then
          write(stringa,'(i11)') lnbc(inpstr)
          write(*,'(4a)')'Error in input file: ',
     &     'ASCII delimitor selection is of length ',
     &     stringa(fnbc(stringa):11),' instead of 3'
          errnum=-2
        elseif((inpstr(1:1).ne.'"').or.(inpstr(3:3).ne.'"')) then
          write(*,'(2a)')'Error in input file: ',
     &     'ASCII delimitor must be enclosed in double-quotes'
          errnum=-2
        else
          delimit=inpstr(2:2)
        endif
      endif
c
c  Read number of min-max pairs for output ASCII file.  A value of zero
c  means that output is a straight conversion of input.
c
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
        write(*,'(a)')'Error in input file: no min-max specified'
        errnum=-2
      elseif(errnum.eq.0) then
        read(inpstr,*,iostat=inpstat) minmax
        if(inpstat.ne.0) then
          write(*,'(a)')
     &     'Error in input file: format error in min-max selection'
          errnum=-2
        elseif(minmax.lt.0) then
          write(stringa,'(i11)') minmax
          write(*,'(3a)')'Error in input file: ',
     &     'invalid minmax specification of ',
     &     stringa(fnbc(stringa):11)
          errnum=-2
        endif
      endif
c
c  Read suntracker total intensity limits (STD and AVG).
c
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
        write(*,'(a)')
     &   'Error in input file: no suntracker intensity limits'
        errnum=-2
      elseif(errnum.eq.0) then
        read(inpstr,*,iostat=inpstat) stlimstd,stlimavg
        if(inpstat.ne.0) then
          write(*,'(2a)')
     &     'Error in input file: format error in ',
     &     'suntracker intensity limits'
          errnum=-2
        elseif((stlimstd.lt.0.0).or.(stlimavg.lt.0.0)) then
          write(*,'(2a)')'Error in input file: ',
     &     'suntracker intensity limits must be positive'
          errnum=-2
        endif
      endif

c
c  Read minimum value of the ZPD peaks.
c
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
        write(*,'(a)')
     &   'Error in input file: no minimum values for the ZPD peaks'
        errnum=-2
      elseif(errnum.eq.0) then
        read(inpstr,*,iostat=inpstat)(zpdlim(indexa),indexa=1,dtc2)
        if(inpstat.ne.0) then
          write(*,'(2a)')
     &     'Error in input file: format error in ',
     &     'minimum values for the ZPD peaks'
          errnum=-2
        else
          do indexa=1,dtc2
            if((errnum.eq.0).and.(zpdlim(indexa).lt.0.0)) then
              write(*,'(2a)')'Error in input file: ',
     &         'ZPD minimum values must be positive'
              errnum=-2
            endif
          enddo
        endif
      endif
c
c  Read minimum value of the cross-correlation between runs.
c
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
        write(*,'(a)')
     &   'Error in input file: no values for minimum cross-correlation'
        errnum=-2
      elseif(errnum.eq.0) then
        read(inpstr,*,iostat=inpstat)(xcorlim(indexa),indexa=1,dtc2)
        if(inpstat.ne.0) then
          write(*,'(2a)')
     &     'Error in input file: format error in ',
     &     'values for minimum cross-correlation'
          errnum=-2
        else
          do indexa=1,dtc2
            if((errnum.eq.0).and.(xcorlim(indexa).lt.0.0)) then
              write(*,'(2a)')'Error in input file: ',
     &         'cross-correlation limits must be positive'
              errnum=-2
            endif
          enddo
        endif
      endif

c
c  Read the time correction, in floating-point hours.
c
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
        write(*,'(a)')
     &   'Error in input file: no time correction'
        errnum=-2
      elseif(errnum.eq.0) then
        read(inpstr,*,iostat=inpstat) timecorr
        if(inpstat.ne.0) then
          write(*,'(a)')
     &     'Error in input file: format error in time correction'
          errnum=-2
        else
          if((errnum.eq.0).and.
     &     ((timecorr.lt.(-24.0)).or.(timecorr.gt.24.0))) then
              write(*,'(2a)')'Error in input file: ',
     &         'time correction must be between -24 and +24'
            errnum=-2
          endif
        endif
      endif
c
c  Read maximum number of scans to analyze.
c
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
        write(*,'(a)')'Error in input file: no maximum number of scans'
        errnum=-2
      elseif(errnum.eq.0) then
        read(inpstr,*,iostat=inpstat) scanlim
        if(inpstat.ne.0) then
          write(*,'(2a)')
     &     'Error in input file: format error in ',
     &     'maximum number of scans'
          errnum=-2
        elseif(scanlim.lt.0) then
          write(stringa,'(i11)') scanlim
          write(*,'(3a)')'Error in input file: ',
     &     'invalid maximum number of scans of ',
     &     stringa(fnbc(stringa):11)
          errnum=-2
        endif
      endif
c
c  Read maximum log-base-2 of the FFT size.
c
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
        write(*,'(a)')
     &   'Error in input file: no values for maximum FFT size'
        errnum=-2
      elseif(errnum.eq.0) then
        read(inpstr,*,iostat=inpstat)(fftlim(indexa),indexa=1,dtc2) 
        if(inpstat.ne.0) then
          write(*,'(2a)')
     &     'Error in input file: format error in ',
     &     'values of maximum FFT size'
          errnum=-2
        else
          do indexa=1,dtc2
            if((errnum.eq.0).and.(fftlim(indexa).lt.0)) then
              write(*,'(2a)')'Error in input file: ',
     &         'log-base-2 FFT size limits must be positive'
              errnum=-2
            endif
          enddo
        endif
      endif
c
c  Read SIV_correction frequencies
c
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
        write(*,'(a)')
     &   'Error in input file: no values for SIVCF'
        errnum=-2
      elseif(errnum.eq.0) then
        read(inpstr,*,iostat=inpstat)(sivcfreq(indexa),indexa=1,dtc2) 
        if(inpstat.ne.0) then
          write(*,'(2a)')
     &     'Error in input file: format error in ',
     &     'values of SIV correction frequencies'
          errnum=-2
        else
          do indexa=1,dtc2
            if((errnum.eq.0).and.(sivcfreq(indexa).lt.0)) then
              write(*,'(2a)')'Error in input file: ',
     &         'values of SIV correction frequencies must be positive'
              errnum=-2
            endif
          enddo
        endif
      endif
c
c  Read LSE type
c
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)')
     &   'Error in input file: no laser sampling error type'
         errnum=-2
      elseif(errnum.eq.0) then
         read(inpstr,*,iostat=inpstat)(lsemode(indexa),indexa=1,dtc2)
         if(inpstat.ne.0) then
            write(*,'(2a)')
     &      'Error in input file: format error in ',
     &      'laser sampling error type'
            errnum=-2
         else
           do indexa=1,dtc2
            if((errnum.eq.0).and.
     &      ((lsemode(indexa).lt.0).or.(lsemode(indexa).gt.4))) then
             write(stringa,'(i11)') lsemode(indexa)
             write(*,'(3a)')'Error in input file: ',
     &       'invalid laser sampling error type of ',
     &       stringa(fnbc(stringa):11)
             errnum=-2
            endif
           enddo
         endif
      endif

c
c  Read lengths of the PCO
c
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
        write(*,'(a)') 'Error in input file: no PCO thresholds'
        errnum=-2
      elseif(errnum.eq.0) then
        read(inpstr,*,iostat=inpstat)(pco_len(indexa),indexa=1,dtc2) 
        if(inpstat.ne.0) then
          write(*,'(2a)')
     &     'Error in input file: format error in PCO lengths'
          errnum=-2
        else
          do indexa=1,dtc2
            if((errnum.eq.0).and.(pco_len(indexa).lt.0)) then
              write(*,'(2a)')'Error in input file: ',
     &         'values of PCO lengths must be positive'
              errnum=-2
            endif
          enddo
        endif
      endif
c
c  Read spectral magnitude thresholds for the PCO
c
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
        write(*,'(a)') 'Error in input file: no PCO thresholds'
        errnum=-2
      elseif(errnum.eq.0) then
        read(inpstr,*,iostat=inpstat)(pco_thresh(indexa),indexa=1,dtc2) 
        if(inpstat.ne.0) then
          write(*,'(2a)')
     &     'Error in input file: format error in PCO thresholds'
          errnum=-2
        else
          do indexa=1,dtc2
            if((errnum.eq.0).and.(pco_thresh(indexa).lt.0)) then
              write(*,'(2a)')'Error in input file: ',
     &         'values of PCO thresholds must be positive'
              errnum=-2
            endif
          enddo
        endif
      endif
c
c  Read maximum level of processing.
c
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
          write(stringa,'(i11)') proclim
          write(*,'(3a)')'Error in input file: ',
     &     'invalid maximum level of processing of ',
     &     stringa(fnbc(stringa):11)
          errnum=-2
        endif
      endif
c
c  Read verbosity level.
c
      call read_input_line(luni,errnum,inpstat,inpstr)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
        write(*,'(a)')'Error in input file: no verbosity specified'
        errnum=-2
      elseif(errnum.eq.0) then
        read(inpstr,*,iostat=inpstat) verbose
        if(inpstat.ne.0) then
          write(*,'(2a)')
     &     'Error in input file: format error in ',
     &     'verbosity specification'
          errnum=-2
        elseif((verbose.lt.1).or.(verbose.gt.5)) then
          write(stringa,'(i11)') verbose
          write(*,'(3a)')'Error in input file: ',
     &     'invalid verbosity specification of ',
     &     stringa(fnbc(stringa):11)
          errnum=-2
        endif
      endif
c
      return
      end
