      program opus_i2s
c
c  Program to convert OPUS interferograms to spectra.
c
c  External parameters and a list of scan sets are obtained from
c  the input file "opus-i2s.in" which is heavily commented.
c
      implicit none
c
c  We first define a version string and a program release date.
c  Please remember to update those when you modify this program.
c
      character
     & verstr*(*)  ! Version string displayed at program start

      real*8
     & progver     ! Program version (date) to be stored in OPUS header

      parameter (verstr=
     & ' opus-i2s     version 2.71      2014-03-02  GCT')
      parameter (progver=20140302.d0)
c
c  The parameters below control how much memory is used by the program.
c
c  The dominant term is MIP which is both the maximum number of
c  interferogram points and the largest FFT size that the program
c  can compute.  Because the program currently duplicates the long
c  side of the interferogram just before the FFT, the actual maximum
c  number of interferogram points is MIP/2 plus the number of points
c  on the short side of ZPD.  For the value of 2**24, the maximum
c  number of points on the long side of ZPD is 8388608.  Given that
c  two points are acquired per laser fringe and that the laser
c  wavenumber is 15798 cm-1, the maximum path difference currently
c  supported is: 8388608 / (2 * 15798) = 265 cm.
c
c FIXME: do we even need MNS and MSL ??
CCc  Parameter MNS is the maximum number of scans in a set.  In the
CCc  IFS125, this number is controlled by NSS (Number of Sample Scans).
CCc  Routinely, NSS will be either 1, 2, or 4.  However the vectors
CCc  controlled by MNS consume little memory, so we can base this choice
CCc  on a continuously scanning instrument (hoping for a firmware upgrade
CCc  that would allow interruptible scanning).  With the current value
CCc  of MNS=100, we could scan for roughly 2.5 hours without interruptions.
c
CCc  In deriving MSL, the maximum number of data slices in a set of scans,
CCc  we have used the maximum total path of (6 + 148) cm and the fact that
CCc  data slices contain about 190000 points.  So the maximum number of
CCc  data slices per run is: (6 + 148) * 2 * 15798 / 190000 = 26.
c
c  Parameter MCH is fixed at 2 because the IFS125 supports only two
c  detectors (master+slave) in its acquisition system.
c
c  Parameter MDTC is the maximum number of detectors that can be handled
c  in a single program execution.  Some information from the input file
c  is detector specific and must be held in appropriately-sized vectors.
c
c  Parameters MI4 and MR8 are the maximum sizes of vectors used to gather
c  all the OPUS header items.
c
c FIXME: this comes from the input file in this pgm
CCc  Parameter MIF is the maximum number of items read from the information
CCc  files produced by the real-time recording system.
c
      integer*4
     & mip,        ! Maximum number of input points or FFT size
     & mns,        ! Maximum number of interferograms per OPUS file
     & msl,        ! Maximum number of interferogram sections per OPUS file
     & mch,        ! Maximum number of data channels
     & mdtc,       ! Maximum number of detectors
     & mi4,        ! Maximum number of I*4 items in file header
     & mr8,        ! Maximum number of R*8 items in file header
     & mif         ! Maximum number of info channels
      parameter (mip=2**24)
      parameter (mns=1)
      parameter (msl=1)
      parameter (mch=2)
      parameter (mdtc=5)
      parameter (mi4=40)
      parameter (mr8=40)
      parameter (mif=40)
c
c  The remaining data declarations are split into two groups.  First here
c  are the variables that the program uses for its internal functions or
c  to store the parameters read from the input file.
c
      integer*4 nfft,
     & pco_len(mdtc),! Length of phase correction operator
     & errnum,      ! Error code (0=ok, <0=fatal, >0=recoverable)
     & inpstat,     ! Value of IOSTAT from parameter input file read
     & lun,         ! Logical Unit Number for file I/O
     & igrmode,     ! Interferogram saving mode (none/raw/deglitched)
     & lsemode(mdtc), ! Laser sampling error type (none/slave/master/Hase/other)
     & phmode,      ! Phase saving mode (none/yes)
     & chan1,       ! Starting channel number to process
     & chan2,       ! Ending channel number to process
     & dchan,
     & outfmt,      ! Format selection for output file
     & minmax,      ! Count of min-max pairs in ASCII output file (0=no min-max)
     & fftlim(mdtc),! Maximum log-base-2 of the FFT size for each detector
     & proclim,     ! Maximum processing stage performed by the program
     & verbose,     ! Level of verbosity for displayed messages
     & catyear,     ! Year from catalog
     & catmonth,    ! Month from catalog
     & catday,      ! Day from catalog
     & runno,       ! Run number increasing throughout the day
     & memrun,      ! Memory of run number for channels > 1
     & scanno,      ! Loop index for separation of FWD and REV
     & scanlim,     ! Not used yet (argument from parse_input_top)
     & scanllt,     ! Loop limit for separation of FWD and REV
     & scancnt,     ! Number of scans in one scanner direction
     & channel,     ! Channel number (1=InGaAs=slave, 2=Si=master)
     & counter,     ! Point count, data destination into 'buf_igram'
     & memcnt,      ! Memory of point count for reverse direction
     & indexa,      ! General loop index
     & lnbc,        ! Integer function Last Non-Blank Character in string
     & fnbc,        ! Integer function First Non-Blank Character in string
     & fbc,         ! Integer function First Blank Character in string
     & lf,lb,       ! Indices inro the input string
     & nspectra,    ! number of spectra in catalogue
     & nbadspectra  ! number of spectra excluded by SIS/SIA criterion

      real*4 
     & frsp,        ! The fraction of the spectral domain without optical energy
     & pinv,
     & timecorr,    ! Time correction, added to instrument time --> UT
     & stlimstd,    ! Suntracker rel std dev limit
     & stlimavg,    ! Suntracker avg intensity limit
     & zpdlim(mdtc),! Minimum value of the ZPD peaks
     & xcorlim(mdtc)! Minimum value of the ZPD peaks

      real*8
     & frzpda,      ! Fractional size of ZPD artifact (dip) in smoothed igram
     & sivcfreq(mdtc),  ! SIV-correction frequencies (cm-1)
     & pco_thresh(mdtc) ! SIV-correction frequencies (cm-1)

      character
     & infile*100,  ! Name of program input file
     & inpath*99,   ! Directory path to input OPUS files
     & outpath*99,  ! Directory path for output spectrum files
     & igrmpath*99, ! Directory path for output interferogram files
     & phpath*99,   ! Directory path for output phase files
     & flimit*99,   ! Name of file containing the frequency limits
     & pattern*99,  ! Pattern for CIT file-naming convention
     & srcindic*26, ! List of source indicators for output files
     & dtcigrm*26,  ! List of detector indicators for interferogram files
     & dtcspec*26,  ! List of detector indicators for spectrum files
     & delimit*1,   ! Field delimiter for ASCII output file
     & opusname*99, ! Name of OPUS file containing interferogram
     & inpstr*999,  ! String used to read entries from input file
     & filename*199, ! Full file name for output files
     & phasepath*99, ! Full file name for output files
     & stringa*11,  ! String used to format integer display
     & cc*40

      logical*4
     & filexist     ! Keeps track of file existence

c      parameter (infile='opus-i2s.in')
      parameter (lun=20)
c
c  The second group of data declarations holds the values of the OPUS data
c  and header items.  These are read (or derived) from the interferogram
c  and some are stored in the output interferogram/spectrum.  A few
c  items (e.g.: nss, tpx) have meanings specific to this program, but most
c  header items are generic: those are stored in two vectors to allow for
c  easy expansion without constantly changing the subroutine argument lists.
c  The indices of the header items kept in those vectors are defined in the
c  following include file:
c
      include '../opus-comn/header_indices.inc'

      integer*4
     & nptvec(msl), ! Number of PoinTs, from OPUS header
     & bpdata(msl,mch),! Byte pointers into the data blocks of interferograms
     & nss,         ! Number of Sample Scans in the input interferograms
     & nsubs,
     & tpx,         ! Total Points in X (along OPD) in each scan direction
     & pinl,        ! Peak interferogram location
     & sivcflag,    ! SIV Correction flag
c     & nlimit,
     & istat,
     & izpd,        ! Point index of zero path difference
     & i4head(mi4)  ! Vector to hold the I*4 header items

      real*4 
     & buf_igram(mip),  ! Data storage and FFT buf_igram
     & smoo_igram(mip), ! Storage for smoothed igram
     & revbuf(mip/2),   ! Temporary storage for reverse run from OPUS file
     & shbar,           ! Calculated laser sampling error (LSE)
     & sherr            ! Calculated laser sampling error uncertainty (LSU)

      real*8
     & dclevel,     ! The DC interferogram signal level at ZPD
     & fvsi_calc,   ! Calculated FVSI from the smoothed interferogram
     & zpa,         ! ZPD interferogram amplitude (phase-corrected)
     & r8head(mr8), ! Vector to hold the R*8 header items
     & timvec(mns), ! Time vector contains one entry for each scan
     & minfold,     ! Wavenumber of the start of this spectral fold
     & maxfold,     ! Wavenumber of the start of the next spectral fold
     & infovec(mif),! Information produced by real-time algorithm
     & Tstart,      ! The time for the first point of the first fwd or rev scan
     & Tdur,        ! DUR parameter as defined in original opus file
     & Tscan,       ! actual time for 1 scan, calculated from NPT and VEL
     & Tturn        ! Turnaround time between fwd-rev scans
c
      character*40
     & DTCstr,      ! Character variable to hold DTC description
     & INSstr       ! Character variable to hold INS description

c  Program execution starts here: display version string.
c
      write(*,'(a)')verstr
c
c  Initialize variables.
c
      errnum=0
      inpstat=0
      i4head(i_bfw)=0    ! We never save bad runs: BFW=BBW=0
      i4head(i_bbw)=0 
      sivcflag=0
      dclevel=0.0d0
      fvsi_calc=0.0d0
      shbar=0.0
      sherr=0.0
      infovec=0.d0
c      nlimit=3600
c
c  Check that input file exists, if so open it.
c
      call getarg(1,infile)
      istat=lnbc(infile)
      if(istat.le.0) then
         write(*,*)'Missing input filename on command line '
         write(*,*)'Using the default -- opus-i2s.in '
         infile='opus-i2s.in'
      endif

      inquire(file=infile,exist=filexist,iostat=inpstat)
      if(inpstat.ne.0) then
        errnum=-1
        write(*,'(2a)')'Error: inquire failed on input file ',infile
      elseif(filexist) then
        open(unit=lun,file=infile,status='old',iostat=inpstat)
        if(inpstat.ne.0) then
          errnum=-1
          write(*,'(2a)')'Error: open failed on input file ',infile
        endif
      else
        errnum=-1
        write(*,'(2a)')'Error: cannot find the input file ',infile
      endif
c
c  Parse the top section of the input file containing general parameters.
c
      call parse_input_top(lun,mdtc,errnum,inpstat,inpath,outpath,
     & igrmode,igrmpath,phmode,phpath,chan1,chan2,flimit,pattern,
     & srcindic,dtcigrm,dtcspec,outfmt,delimit,minmax,stlimstd,
     & stlimavg,zpdlim,xcorlim,timecorr,scanlim,fftlim,sivcfreq,lsemode,
     & pco_len,pco_thresh,proclim,verbose)

c
c  Main program loop, over the scan sets listed in the input file.
c
        nspectra=0
        nbadspectra=0
        do while(errnum.eq.0)
        call read_input_line(lun,errnum,inpstat,inpstr)
          if(inpstat.ne.0) exit
          nspectra=nspectra+1
c
c  Reset "Number of Sample Scans" in case there is an error detected
c  before it can be read from the interferogram headers.
c
        nss=0
        scanllt = 1
        scancnt = nss
c
c  Parse one line of the input file.
c
        call substr(inpstr,cc,1,nsubs)
        lb=fnbc(inpstr)
        lf=fbc(inpstr(lb:))+lb
        opusname=inpstr(lb:lf-1)
c        write(*,*)nsubs,lb,lf,opusname
        if(nsubs.eq.16) then
           read(inpstr(lf:),*,iostat=inpstat)
     &      catyear,catmonth,catday,runno,
     &      infovec(7),infovec(8),infovec(9),           ! LAT, LON, ALT
     &      infovec(4),infovec(3),infovec(1),           ! TSC, PIM, HUM
     &      infovec(15),infovec(19),infovec(16),        ! TOU, POU, HOU
     &      infovec(28), infovec(29)                    ! SIA, fvsi (fractional)
        elseif(nsubs.eq.18) then
           read(inpstr(lf:),*,iostat=inpstat)
     &      catyear,catmonth,catday,runno,
     &      infovec(7),infovec(8),infovec(9),           ! LAT, LON, ALT
     &      infovec(4),infovec(3),infovec(1),           ! TSC, PIM, HUM
     &      infovec(15),infovec(19),infovec(16),        ! TOU, POU, HOU
     &      infovec(28), infovec(29),                   ! SIA, fvsi (fractional)
     &      infovec(10), infovec(13)                    ! WSPD, WDIR  added DG 110220
        else
            stop 'unrecognized infor file format'
        endif

        infovec(29)=infovec(29)*infovec(28)          ! convert fvsi to absolute(SIS)
        infovec(2)=infovec(4)                        ! TLP
        if(inpstat.eq.-1)then
           write(*,'(/,a,/,a)')
     &     'End of string reading next input line:
     & some parameters may be missing', lnbc(inpstr)
          errnum=0
          inpstat=0
        elseif(inpstat.ne.0) then
          write(*,'(a)')
     &     'Error in input file: format error in catalog'
          errnum=1
          inpstat=0
        endif
c
        if(stlimavg.ne.0.0.and.stlimstd.ne.0)then
          if(infovec(28).lt.stlimavg.or.
     &      infovec(29).gt.stlimstd*infovec(28))then
            write(*,'(a,a,a,f6.1,a,f6.2)') opusname(:lnbc(opusname)),
     &      ': Tracker intensity too low or stdev too ',
     &      'high: SIA = ', infovec(28), 
     &      ', SIS = ',infovec(29)
              nbadspectra=nbadspectra+1
             cycle
          endif
          endif
c
c  Obtain all IFS12x parameters by reading the header of the OPUS file.
c
        if((errnum.eq.0).and.(proclim.ge.1)) then
          call get_runopus_parameters(inpath,opusname,
     &     runno,verbose,mns,msl,mip,mch,mi4,mr8,chan1,chan2,
     &     errnum,nptvec,bpdata,nss,tpx,timvec,i4head,r8head,
     &     DTCstr,INSstr)
        endif
c
c  Setup for FWD-only and FWD+REV runs.
c
        Tscan=0.0d0
        Tturn=0.0d0
        if((errnum.eq.0).and.(proclim.ge.1)) then
          if((i4head(i_aqmcode).eq.aqm_sf).or.
     &       (i4head(i_aqmcode).eq.aqm_sn)) then
            scanllt = 1
            scancnt = nss
          elseif(i4head(i_aqmcode).eq.aqm_sd) then
            if(nss.gt.1) then
              scanllt = 2
              scancnt = nss/2
c     DG090401:  calculate scan time and turnaround time 
              Tdur = r8head(i_dur)        !Total duration time first-last point
              Tscan = dble(tpx/2)/2.0D0/r8head(i_laserate)
     &               /dble(i4head(i_SSM))*dble(i4head(i_SSP))
              if(nss.gt.1) then
                  Tturn = (Tdur - nss*Tscan)/dble(nss-1)
              else
                  Tturn=0
              endif
              r8head(i_dur) = Tdur - Tscan - Tturn
c     DG090401 end
              if(mod(runno,2).ne.1) then
                errnum=2
                write(*,'(a)')
     &          'Error: run number for output FWD file must be odd'
              endif
              if(mod(tpx,2).ne.0) then
                errnum=2
                write(stringa,'(i11)') tpx
                write(*,'(4a)')'Error: ',opusname(1:lnbc(opusname)),
     &           ' has odd TPX of ',stringa(fnbc(stringa):11)
              else
                tpx=tpx/2
              endif
            else
              scanllt = 1
              scancnt = nss
            endif
          else
            errnum=99
            write(*,'(a)')
     &       'Unimplemented: AQM is different from SF, SN, or SD'
          endif
        endif
c
c  Setup the OPUS header items that will remain unchanged throughout
c  this scan set.
c
        if((errnum.eq.0).and.(proclim.ge.1)) then
c
c  Apply site-specific time correction.
c
          timvec(1)=timvec(1)+(dble(timecorr)*3600.d0)
c
c  FIXME: the zero-filling factor could be set to values greater than 2
c  if we allowed the FFT size to grow beyond the nearest power of two.
c
          i4head(i_zff)=2
c
c  Compute the folding limits for this spectral domain.  We dont use the
c  information from the interferogram headers because it is not accurate
c  enough.  We did make sure that LWN is correct, so we use that instead.
c
          minfold=0.d0
          maxfold=minfold+((r8head(i_lwn)*
     &     dble(i4head(i_ssm)))/dble(i4head(i_ssp)))
        endif
c
c  Now we process each channel for the FWD then REV runs.
c
        do scanno=1,scanllt
        do channel=chan2,chan1,-1
c
c  Read interferogram(s)   (combined FWD+REV if AQM=SD and nss>1)
          if((errnum.eq.0).and.(proclim.ge.3)) then
            call get_runopus_data(inpath,opusname,verbose,
     &       msl,mip,mch,nptvec,bpdata,channel,counter,buf_igram)
            nfft=1
            indexa=0
           do while((nfft.lt.counter/2).and.(indexa.lt.fftlim(channel)))
              nfft=nfft*2
              indexa=indexa+1
            enddo

c
c  If file was recorded in SD mode and contains more than one scan,
c  divide counter by two to separate FWD from REV runs.  Also save
c  REV data into temporary buf_igram.
c
c FIXME: check that TPX and counter match
            if((errnum.eq.0).and.(scanllt.eq.2)) then
              counter=counter/2
c              do indexa=1,counter
c                revbuf(indexa)=buf_igram(counter+indexa)
c              enddo
              do indexa=counter,1,-1
                 buf_igram(nfft+indexa)=buf_igram(counter+indexa)
              end do
            endif
          endif
c
c  Now loop over the run directions (FWD-only or FWD+REV)
          if((errnum.eq.0).and.(proclim.ge.3)) then
c            memrun=runno
c            memcnt=counter
c
c  Set run direction variable based on parity of run number (odd=fwd).
c  For reverse runs, get data back from temporary buf_igram.
c  DG090226 set TIM as start time for fwd and rev scans 
c

c     DG090522 The next line crashes Lauder runs with only 
c     FWD single scans and consecutive run numbers
c              if(mod(runno,2).ne.0) then  !FWD scan
c     Replace with:
              if(scanllt.eq.1) then !single scans
                if(i4head(i_gbw).eq.1) then ! SLICE-I2S REV scans
                  i4head(i_gfw)=0
c SLICE-I2S generated reverse scan: want buf_igram, not revbuf
                  Tstart = timvec(1) + Tscan + Tturn
                 if(i4head(i_prl).lt.nptvec(1)/3) then
c This implies that the backward interferogram has already been
c reversed, looks like a forward interferogram and should be treated as
c such. The divisor (3) is somewhat arbitrary. It should be >=2.
                   i4head(i_gbw)=0
                   i4head(i_gfw)=1
                   write(*,*)'Treating as FWD scan.'
                 endif
                else ! FWD scans
                  i4head(i_gfw)=scancnt
                  Tstart=timvec(1)
                endif
              else
                if(mod(runno,2).ne.0) then  !FWD scan
                i4head(i_gfw)=scancnt
                Tstart=timvec(1)
              else                        !REV scan
                i4head(i_gfw)=0
c                  do indexa=1,counter
c                    buf_igram(indexa)=revbuf(indexa)
c                  enddo
c              do indexa=1,counter
c                revbuf(indexa)=buf_igram(counter+indexa)
c              enddo
                  do indexa=counter,1,-1
                    buf_igram(indexa)=buf_igram(counter+indexa)
                  end do
                Tstart = timvec(1) + Tscan + Tturn
              endif
              endif
              i4head(i_gbw)=scancnt-i4head(i_gfw)
c
c  If enabled save separated interferograms in the requested format.
c
              if((errnum.eq.0).and.(proclim.ge.3)) then
                if(igrmode.eq.1) then
                  call build_cit_name(igrmpath,pattern,dtcigrm,srcindic,
     &             catyear,catmonth,catday,channel,i4head(i_srccode),
     &             runno,filename)
                  call save_to_file(1,flimit,filename,channel,outfmt,
     &             verbose,progver,mip,mif,mi4,mr8,
     &             buf_igram(1+nfft*(scanno-1)),counter,
     &             minfold,maxfold,delimit,minmax,tpx,Tstart,
     &             i4head,r8head,DTCstr,INSstr,sivcfreq(channel),
     &             pco_len(channel),pco_thresh(channel),izpd,
     &             sivcflag,dclevel,fvsi_calc,zpa,
     &             shbar,sherr,lsemode(channel),
     &             infovec,tla_ext_min,errnum)
                endif
              endif
c
              frsp=(sivcfreq(channel)-minfold)/(maxfold-minfold)
c              write(*,*)'Calling SIV correction'
              call siv_correction(counter,buf_igram(1+nfft*(scanno-1)),
     &        smoo_igram,
     &        frsp,pinl,pinv,sivcflag,dclevel,fvsi_calc,frzpda)
c              write(*,*)'Calling SIV correction'

cc  Dont let ZPD be too close to ends of interferograms
cc  Such cases tend to be garbage scans anyway and can cause
cc  very small OPD values in the runlog which cause GFIT problems
c            if (pinl.lt.nlimit) pinl=nlimit
c            if (pinl.gt.counter-nlimit) pinl=counter-nlimit

c
c  Perform Fourier transform and save data in the requested format.
c
              if((errnum.eq.0).and.(proclim.ge.4)) then
                if(phmode.gt.0) then
                  call build_cit_name(phpath,pattern,dtcspec,srcindic,
     &             catyear,catmonth,catday,channel,i4head(i_srccode),
     &             runno,phasepath)
                   phasepath=phasepath(:lnbc(phasepath))//'.phs'
                else
                  phasepath=''
                endif

                call i2s_lite_siv(mip,fftlim(channel),verbose,
     &          i4head(i_gfw),pinl,
     &          channel,pco_len(channel),pco_thresh(channel),
     &          phasepath,buf_igram(1+nfft*(scanno-1)),counter,izpd,zpa,
     &          shbar,sherr,lsemode(channel))

                call build_cit_name(outpath,pattern,dtcspec,srcindic,
     &           catyear,catmonth,catday,channel,i4head(i_srccode),
     &           runno,filename)

                call save_to_file(2,flimit,filename,channel,outfmt,
     &          verbose,progver,mip,mif,mi4,mr8,
     &          buf_igram(1+nfft*(scanno-1)),counter,
     &           minfold,maxfold,delimit,minmax,tpx,Tstart,i4head,
     &           r8head,DTCstr,INSstr,sivcfreq(channel),
     &           pco_len(channel),pco_thresh(channel),izpd,sivcflag,
     &           dclevel,fvsi_calc,zpa,shbar,sherr,lsemode(channel),
     &           infovec,tla_ext_min,errnum)
              endif
c              counter=memcnt
          endif
        enddo      ! channel=chan1,chan2
              runno=runno+1
c            runno=memrun
            enddo    ! scanno=1,scanllt
c
c  Positive error codes indicate transient errors, specific to one set
c  of scans.  The program can proceed to the next set of scans.
c
        if(errnum.gt.0) then
          errnum=0
        endif
c
      enddo      ! while (errnum.eq.0))
c
c  Close the input file
c
      close(unit=lun,iostat=inpstat)
      if(inpstat.ne.0) then
        write(*,'(2a)')'Error: close failed on input file ',infile
      endif

      print *
      print '(2(I6,a))', nspectra, ' opus files in catalogue, ',
     &                   nbadspectra, ' excluded by SIA/SIS'
c
c  Terminate.
c
      stop
      end
c
c ===================================================================
      subroutine get_runopus_parameters(inpath,opusname,
     & runno,verbose,mns,msl,mip,mch,mi4,mr8,chan1,chan2,
     & errnum,nptvec,bpdata,nss,tpx,timvec,i4head,r8head,DTCstr,INSstr)
c
c  Input:
c    inpath    C*(*)  Directory path to input OPUS file
c    opusname  C*(*)  Name of input OPUS file
c    runno     I*4    Run number increasing throughout the day
c    verbose   I*4    Level of verbosity for displayed messages
c    mns       I*4    Maximum number of interferograms per scan set (max NSS)
c    msl       I*4    Maximum number of interferogram slices per scan set
c    mip       I*4    Maximum number of input points
c    mch       I*4    Maximum number of data channels
c    mi4       I*4    Maximum number of I*4 items in file header
c    mr8       I*4    Maximum number of R*8 items in file header
c    chan1     I*4    Starting channel number to process
c    chan2     I*4    Ending channel number to process
c
c  Input/output:
c    errnum    I*4    Error code (0=ok, <0=fatal, >0=recoverable)
c
c  Output:
c    nptvec(msl)I*4   Number of PoinTs, from OPUS header
c    bpdata(msl,mch)I*4 Byte pointers into the data blocks of input file
c    nss       I*4    Number of Sample Scans in this set
c    tpx       I*4    Number of points in one FWD or REV scan
c    timvec(mns)R*8   Time vector contains one entry for each scan
c    i4head(mi4)I*4   Vector to hold the I*4 header items
c    r8head(mr8)R*8   Vector to hold the R*8 header items
c    DTCstr     C*40  Variable to hold detector description
c
      implicit none

      integer*4
     & runno,      ! Subroutine input argument (see above)
     & verbose,    ! Subroutine input argument (see above)
     & mns,        ! Subroutine input argument (see above)
     & msl,        ! Subroutine input argument (see above)
     & mip,        ! Subroutine input argument (see above)
     & mch,        ! Subroutine input argument (see above)
     & mi4,        ! Subroutine input argument (see above)
     & mr8,        ! Subroutine input argument (see above)
     & chan1,      ! Subroutine input argument (see above)
     & chan2,      ! Subroutine input argument (see above)
     & errnum,     ! Subroutine input/output argument (see above)
     & nptvec(msl),! Subroutine output argument (see above)
     & bpdata(msl,mch),! Subroutine output argument (see above)
     & nss,        ! Subroutine output argument (see above)
     & tpx,        ! Subroutine output argument (see above)
     & i4head(mi4),! Subroutine output argument (see above)
     & filecnt,    ! Junk variable to keep get_opusigram_params happy
     & lnbc,       ! Integer function Last Non-Blank Character in string
     & fnbc,       ! Integer function First Non-Blank Character in string
     & scancnt,    ! Count of scans
     & counter,    ! Count of data points
     & errmem      ! Value of errnum at start of scan splitting

      real*8
     & timvec(mns),! Subroutine output argument (see above)
     & r8head(mr8),! Subroutine output argument (see above)
     & timsli(2600)! Time vector contains one entry for each slice
c FIXME: the 2600 should really be MSL, but Sun Fortran hates that

      character
     & inpath*(*), ! Subroutine input argument (see above)
     & opusname*(*),! Subroutine input argument (see above)
     & filename*199,! Full file name for input OPUS file
     & stringa*11, ! String used to format integer display
     & DTCstr*(*), ! Character variable to hold DTC description
     & INSstr*(*)  ! Character variable to hold INS description


      logical*4
     & filexist,   ! Keeps track of file existence
     & inrun       ! Junk

      include '../opus-comn/header_indices.inc'
c
c  Write a blank line to the screen to separate scan sets.
c
      if(verbose.ge.3) then
        write(*,'("")')
      endif
c
c  Check that the OPUS input file exists.
c
c      write(*,*)'inpath=',inpath
c      write(*,*)'opusname=',opusname
      write(filename,'(2a)')inpath(1:lnbc(inpath)),
     &                      opusname(1:lnbc(opusname))
      inquire(file=filename,exist=filexist)
      if(.not.filexist) then
        errnum=1
        write(*,'(3a)')'Error: OPUS file ',
     &   filename(1:lnbc(filename)),' does not exist'
      endif
c
c  Initialize variables.
c
      inrun=filexist
      filecnt=0
      i4head(i_sfmcode)=0  ! Must be initialized: it is used to reject scans.

      if (errnum.eq.0) then
c
c  Read full set of parameters from this OPUS file.
c  Note: 'get_opusigram_params' sets 'inrun' to false at end of scan set.
c
        call get_opusigram_params(filename,verbose,msl,mch,mi4,mr8,
     &   chan1,chan2,errnum,filecnt,inrun,nptvec,bpdata,timsli,nss,tpx,
     &   i4head,r8head,DTCstr,INSstr)

      endif
c
c  If file covers more than one scan, divide TPX by two to separate
c  FWD from REV runs.
c
c     if((errnum.eq.0).and.(nss.gt.1)) then
c       if(mod(tpx,2).ne.0) then
c         errnum=2
c         write(stringa,'(i11)') tpx
c         write(*,'(4a)')'Error: ',filename(1:lnbc(filename)),
c    &     ' has odd TPX of ',stringa(fnbc(stringa):11)
c       else
c         tpx=tpx/2
c       endif
c     endif
c
c  Make sure TPX does not exceed MIP, the maximum number of input points.
c
      if((errnum.eq.0).and.(tpx.gt.mip)) then
        errnum=2
        write(stringa,'(i11)') tpx
        write(*,'(2a)')'Error: increase parameter MIP to ',
     &   stringa(fnbc(stringa):11)
      endif
c
c  Fix laser wavenumber: IFS125-1 instrument returns an incorrect value
c
      if((errnum.eq.0).and.
     & (dabs(r8head(i_lwn)-15798.1d0).lt.(0.01d0))) then
        r8head(i_lwn)=15798.0138d0
        if(verbose.ge.3) then
          write(*,'(a,f11.4)')
     &     'Replacing laser wavenumber of 15798.1 with',r8head(i_lwn)
        endif
      endif
c
c  Initialize scan-splitting variables.
c
      counter=0
      scancnt=0
      errmem=errnum
c
c  FIXME: Unnecessary loop over the slices in this scan set
c         to separate the individual scans.
c
      if(errnum.eq.0) then
        counter=counter+nptvec(1)
        if((counter.eq.tpx).and.(scancnt.eq.mns)) then
          errnum=5
          write(*,'(a)') 'Error: increase parameter MNS'
        elseif(counter.eq.tpx) then
          if(verbose.ge.3) then
            write(stringa,'(i11)') runno+scancnt
            write(*,'(4a)') 'Run ',stringa(fnbc(stringa):11),
     &       ' found'
          endif
          scancnt=scancnt+1
          timvec(scancnt)=timsli(1)
          counter=0
        endif
      endif
c
c  If there was an error, flag undetected runs as invalid.
c
c      do while (scancnt.lt.nss)
c        scancnt=scancnt+1
c        timvec(scancnt)=0.d0
c      enddo
c
c  Restore error state to its value before scan splitting so that runs
c  before the error condition (if any) can be processed.
c
      errnum=errmem

      return
      end
c ===================================================================
      subroutine get_runopus_data(inpath,opusname,verbose,
     & msl,mip,mch,nptvec,bpdata,channel,counter,buf_igram)
c
c  Input:
c    inpath    C*(*)  Directory path to input OPUS files
c    opusname  C*(*)  Name of input OPUS file
c    verbose   I*4    Level of verbosity for displayed messages
c    msl       I*4    Maximum number of interferogram slices per scan set
c    mip       I*4    Maximum number of input points
c    mch       I*4    Maximum number of data channels
c    nptvec(msl)I*4   Number of PoinTs, from OPUS header
c    bpdata(msl,mch)I*4 Byte pointers into the data blocks of slices
c    channel   I*4    Channel number (1=InGaAs=slave, 2=Si=master)
c
c  Output:
c    counter   I*4    Point count, data destination into 'buf_igram'
c    buf_igram(mip)R*4   Data storage and FFT buf_igram
c
      implicit none

      integer*4
     & verbose,    ! Subroutine input argument (see above)
     & msl,        ! Subroutine input argument (see above)
     & mip,        ! Subroutine input argument (see above)
     & mch,        ! Subroutine input argument (see above)
     & nptvec(msl),! Subroutine input argument (see above)
     & bpdata(msl,mch),! Subroutine input argument (see above)
     & channel,    ! Subroutine input argument (see above)
     & counter,    ! Subroutine output argument (see above)
     & iendian,    ! Endianness of computer
     & bigendian,  ! Named constant for big endian detection 
     & luns,       ! Logical Unit Number for OPUS input file
     & reclen,     ! Record length in bytes for single values
     & blklen,     ! Record length in bytes for data blocks
     & blkpnt,     ! Number of data points in a block
     & lnbc,       ! Integer function Last Non-Blank Character in string
     & fnbc,       ! Integer function First Non-Blank Character in string
     & npts,       ! Number of points in OPUS file
     & inpoint,    ! Byte pointer to data values in OPUS file
     & rs1,        ! Starting record number single value reads at start
     & re1,        ! Ending record number single value reads at start
     & rs2,        ! Starting record number block reads
     & re2,        ! Ending record number block reads
     & rs3,        ! Starting record number single value reads at end
     & re3,        ! Ending record number single value reads at end
     & recnum,     ! Number of current data record from OPUS file
     & indexa      ! Index into buf_igram during block read and byte swapping

      parameter (bigendian=1)
      parameter (luns=21)
      parameter (reclen=4)
      parameter (blklen=2**9)
      parameter (blkpnt=blklen/reclen)

      real*4 buf_igram(mip) ! Subroutine output argument (see above)

      character
     & inpath*(*), ! Subroutine input argument (see above)
     & opusname*(*),! Subroutine input argument (see above)
     & filename*199,! Full file name for OPUS file
     & stringa*11, ! String used to format integer display
     & stringb*11  ! String used to format integer display
c
c  Initialize variables.
c
      call getendian(iendian)
      counter=0
 
      inpoint=bpdata(1,channel)
      npts=nptvec(1)

      write(filename,'(2a)')inpath(1:lnbc(inpath)),
     &                      opusname(1:lnbc(opusname))

      if(verbose.ge.4) then
        write(stringa,'(i11)') npts
        write(stringb,'(i11)') inpoint
        write(*,'(5a)') filename(1:lnbc(filename)),
     &   '  npts=',stringa(fnbc(stringa):11),
     &   '  byte=',stringb(fnbc(stringb):11)
      endif
c
c  We pre-compute the limits of three zones for reading the data.  Zone 1
c  is read in one-point increments until we get to a block boundary.
c  Zone 2 is read in full blocks.  Zone 3 goes back to one-point reads
c  to the end of the data section.
c
c  First compute rsX and reX in units of bytes, numbered from 0.
c
      rs1=inpoint
      rs2=(inpoint/blklen)*blklen
      if(rs2.lt.inpoint) rs2=rs2+blklen
      rs3=((inpoint+(npts*reclen))/blklen)*blklen

      re1=rs2-1
      re2=rs3-1
      re3=(inpoint+(npts*reclen))-1 
c
c  Now convert rsX and reX to blocks, numbered from 1.
c
      rs1=1+(rs1/reclen)
      re1=1+(re1/reclen)
      rs2=1+(rs2/blklen)
      re2=1+(re2/blklen)
      rs3=1+(rs3/reclen)
      re3=1+(re3/reclen)
c
c  Convert sample count to a pointer into an array numbered from 1.
c
      counter=counter+1
c
c  Perform file read in three zones (points/blocks/points).
c
      open(luns,file=filename,form='unformatted',status='old',
     & access='direct',recl=reclen) 
      do recnum=rs1,re1
        read(luns,rec=recnum) buf_igram(counter)
        counter=counter+1
      enddo    ! recnum=rs1,re1
      close(luns)

      open(luns,file=filename,form='unformatted',status='old',
     & access='direct',recl=blklen) 
      do recnum=rs2,re2
        read(luns,rec=recnum)
     &   (buf_igram(indexa),indexa=counter,counter+blkpnt-1)
        counter=counter+blkpnt
      enddo    ! recnum=rs2,re2
      close(luns)

      open(luns,file=filename,form='unformatted',status='old',
     & access='direct',recl=reclen) 
      do recnum=rs3,re3
        read(luns,rec=recnum) buf_igram(counter)
        counter=counter+1
      enddo    ! recnum=rs3,re3
      close(luns)
c
c  Convert pointer back to sample count.
c
      counter=counter-1
c
c  On big-endian computer (e.g. SPARC and G4) swap data bytes.
c
      if(iendian.eq.bigendian) then
        do indexa=counter-npts+1,counter
          call r4rev(buf_igram(indexa))
        enddo
      endif

      return
      end
c ===================================================================
c      include 'opus-comn/read_input_line.f'
c      include 'opus-comn/get_opus_xx.f'
c      include 'opus-comn/build_cit_name.f'
c      include 'opus-comn/save_to_file.f'
c      include 'comn/lnbc.f'
c      include 'comn/fnbc.f'
c      include 'comn/getendian.f'
c      include 'comn/xxrev.f'
c      include 'comn/julian.f'
c      include 'opus-comn/get_opusigram_params.f'

c      include 'i2s-lite-siv.f'
c-ftnchek      include '../misc/i2s-lite-stub.f'
