      Program slice_i2s

c  Program to convert OPUS-format interferogram files/slices to spectra.

c  External parameters and a list of scan sets are obtained from
c  the input file "xxxx-i2s.in" which is heavily commented.

c  This program also reads the information files produced by the
c  real-time software, some actions (e.g. scan rejection) are based
c  on that information, and most of it is saved in the OPUS file
c  headers.

      implicit none

      character
     & verstr*44  ! Version string displayed at program start

      real*8
     & progver    ! Program version (date) to be stored in OPUS header

c  We first define a version string and a program release date.
c  Please remember to update those when you modify this program.
      parameter (verstr='i2s     version 3.21       2019-11-13    GCT')
      parameter (progver=20191113.d0)

c  The parameters below control how much memory is used by the
c  program, keeping in mind that memory is allocated statically.

c  The dominant term is MIP which is both the maximum number of
c  interferogram points and the largest FFT size that the program
c  can compute.  Because the program currently duplicates the long
c  side of the interferogram just before the FFT, the actual maximum
c  number of interferogram points is MIP/2 plus the number of points
c  on the short side of ZPD.  For the value of 2**24, the maximum
c  number of points on the long side of ZPD is 8388608.  Given that
c  two points are acquired per laser fringe and that the laser
c  wavenumber is 15798 cm-1, the maximum path difference currently
c  supported is: 8388608 / (2 * 15798) = 265 cm, which exceeds the
c  capability of the interferometer (148 cm).

c  Parameter MNS is the maximum number of scans in a set.  In the
c  IFS125, this number is controlled by NSS (Number of Sample Scans).
c  Routinely, NSS will be either 1, 2, or 4.  However the vectors
c  controlled by MNS consume little memory, so we can base this choice
c  on a continuously scanning instrument (hoping for a firmware upgrade
c  that would allow interruptible scanning).  With the current value
c  of MNS=100, we could scan for roughly 2.5 hours without interruptions.

c  In deriving MSL, the maximum number of data slices in a set of scans,
c  we have used the maximum total path of (6 + 148) cm and the fact that
c  data slices contain about 190000 points.  So the maximum number of
c  data slices per run is: (6 + 148) * 2 * 15798 / 190000 = 26.
c  At 45 cm OPD used by TCCON, there are 9 slices per scan, and 18 slices
c  for a FWD-REV pair.

c  Parameter MCH is fixed at 2 because the IFS125 supports only two
c  detectors (master+slave) in its acquisition system.

c  Parameter MDTC is the maximum number of detectors that can be handled
c  in a single program execution.  Some information from the input file
c  is detector specific and must be held in appropriately-sized vectors.

c  Parameters MI4 and MR8 are the maximum sizes of vectors used to gather
c  all the OPUS header items.

c  Parameter MIF is the maximum number of items read from the information
c  files produced by the real-time recording system.

      integer*4
     & mip,        ! Maximum number of input points or FFT size
     & mns,        ! Maximum number of igrams per scan set
     & msl,        ! Maximum number of igram slices per scan set
     & mch,        ! Maximum number of data channels
     & mdtc,       ! Maximum number of detectors
     & mi4,        ! Maximum number of I*4 items in file header
     & mr8,        ! Maximum number of R*8 items in file header
     & mif         ! Maximum number of info channels

      parameter (mip=2**24)
      parameter (mns=100)
      parameter (msl=mns*26)
      parameter (mch=2)
      parameter (mdtc=5)
      parameter (mi4=40)
      parameter (mr8=40)
      parameter (mif=40)

c  The remaining data declarations are split into two groups.  First here
c  are the variables that the program uses for its internal functions or
c  to store the parameters read from the input file "xxxxx-i2s.in".
      integer*4 inlen,tla_ext,idum,i,
     & nffthr,nside,nburst,margin,
     & pco_len(mdtc), ! Length of phase correction operator
     & errnum,        ! Error code (0=ok, <0=fatal, >0=recoverable)
     & inpstat,       ! Value of IOSTAT from parameter input file read
     & lunr_in,       ! Logical Unit Number for file I/O
     & igrmode,       ! Interferogram saving mode (none/raw/deglitched)
     & lsemode(mdtc), ! Laser sampling error type (none/slave/master/Hase/other)
     & phmode,        ! Phase saving mode (none/yes)
     & chan1,         ! Starting channel number to process
     & chan2,         ! Ending channel number to process
c     & dchan,
     & outfmt,        ! Format selection for output file
     & minmax,        ! Count of min-max pairs in ASCII output file (0=no min-max)
     & mscan,         ! Max number of scans analyzed in a single execution
     & fftlim(mdtc),  ! Max size of phase-corrected igram
     & proclim,       ! Max processing stage performed by the program
     & verbose,       ! Level of verbosity for displayed messages
     & catyear,       ! Year from catalog
     & catmonth,      ! Month from catalog
     & catday,        ! Day from catalog
     & catbatch,      ! Batch (save set) number from catalog
     & catslice,      ! Slice number from catalog, used in slice file name
c     & slicmem,       ! Previous slice number for catalog work-around
c     & memyear,       ! Memory of year to control reset of run number
c     & memmonth,      ! Memory of month to control reset of run number
c     & memday,        ! Memory of day to control reset of run number
     & run_start,     ! Starting run number from input header
     & runno,         ! Run number increasing throughout the day
     & iscan,         ! Scan number within a set of NSS scans
     & scanllt,       ! Loop limit for separation of FWD and REV
     & scancnt,       ! Number of scans analyzed in each direction
     & ichan,         ! Channel number (1=InGaAs=slave, 2=Si=master)
     & nip,           ! Igram Point count, data destination into 'buf_igram'
     & runsta(mns),   ! Run starting slice number (if > 0) or run error (if < 0)
     & runend(mns),   ! Run ending slice number
     & lnbc,          ! Integer function Last Non-Blank Character
     & fnbc,          ! Integer function First Non-Blank Character
     & fbc,           ! Integer function First Blank Character
     & lf,lb,         ! Indices in the input string
     & nil,           ! number of input lines (spectra or slices)
     & nbad,          ! # bad slices/spectra excluded by SIA/SIS criterion
     & ia             ! General loop index
     
      real*4 
     & frsp,          ! Fraction of spectral domain without optical energy
     & pinv,          ! Peak INterferogram Value
     & stlimavg,      ! Limit for suntracker intensity (average)
     & stlimstd,      ! Limit for suntracker intensity (standard deviation)
     & ylimits(mch,2),! Limits of allowed igram values (Min/Max,Master/Slave)
     & xcorlim(mdtc)  ! Minimum value of the ZPD cross-correlation

      real*8
     & frzpda,         ! Fractional size of ZPD artifact (dip) in smoothed igram
     & sivcfreq(mdtc), ! SIV-correction frequencies (cm-1)
     & pco_thresh(mdtc) ! SIV-correction frequencies (cm-1)

      character
     & infile*100,   ! Name of program input file
     & inpath*128,    ! Directory path to OPUS-format input slices/files
     & outpath*128,   ! Directory path for output spectrum files
     & igrmpath*128,  ! Directory path for output interferogram files
     & flimit*128,    ! Name of file containing the frequency limits
     & pattern*128,   ! Pattern for CIT file-naming convention
     & srcindic*26,  ! List of source indicators for output files
     & dtcigrm*26,   ! List of detector indicators for igram files
     & dtcspec*26,   ! List of detector indicators for spectrum files
     & delimit*1,    ! Field delimiter for ASCII output file
     & inpstr*999,   ! String used to read entries from input file
     & phpath*128,   ! Directory path for phase output files
     & phasefull*128,! Full path+file name for phase output files
     & fpsfname*128, ! (First Part of) filename containing OPUS igrams
     & filename*128, ! Full file name for output files
     & cc*64

      logical*4
     & filexist      ! Keeps track of file existence

      parameter (lunr_in=20,
     &            nburst=15)

c  The second group of data declarations holds the values of the OPUS data
c  and header items.  These are read (or derived) from the interferogram
c  slices and some are stored in the output interferogram/spectrum.  A few
c  items (e.g.: nss, tpx) have meanings specific to this program, but most
c  header items are generic: those are stored in two vectors to allow for
c  easy expansion without constantly changing the subroutine argument lists.
c  The indices of the header items kept in those vectors are defined in the
c  following include file:
      include 'header_indices.inc'

      integer*4 iyyy,iyyywas,im,imwas,id,idwas,jd,j0,isi,
     & nptvec(msl),     ! Number of PoinTs in each slice/scan, from OPUS header
     & bpdata(msl,mch), ! Byte pointers into the data blocks of igrams/slices
     & nss,             ! Number of Sample Scans in the input igrams/slices
     & nsubstr,         ! Number of sub-strings: slice-i2s=5; opus-i2s=16 or 18
     & tpx,             ! Total Points in X (along OPD) in each scan direction
     & nlong,           ! # points on long side of phase-corrected igram
     & nshort,          ! # points on short side of phase-corrected igram
     & pinl,            ! Peak interferogram location
     & izpd,            ! Point index of zero path difference
     & sivcflag,        ! SIV correction flag (1 when performed, otherwise 0)
     & istat,           ! Status Flag
     & i4head(mi4)      ! Vector to hold the I*4 header items

      real*4
     & buf_igram(mip),  ! Data storage and FFT buf_igram
     & smoo_igram(mip), ! Storage for smoothed igram
     & timecorr,        ! Time correction (hours) to be added to instrument time
     & shbar,           ! Calculated laser sampling error (LSE)
     & sherr,           ! Calculated laser sampling error uncertainty (LSU)
     & fpilha           ! Fraction of Power in Lower Half of Alias

      real*8
     & dclevel,         ! The DC interferogram signal level at ZPD
     & fvsi_calc,       ! Calculated FVSI from the smoothed interferogram
     & zpa,             ! ZPD interferogram amplitude (phase-corrected)
     & r8head(mr8),     ! Vector to hold the R*8 header items
     & Tstart,
     & timvec(mns),     ! Time vector contains one entry for each scan
     & timsli(msl),     ! Time vector contains one entry for each slice
     & minfold,         ! Wavenumber of the start of this spectral fold
     & maxfold,         ! Wavenumber of the start of the next spectral fold
     & infovec(mif),    ! Information produced by real-time algorithm
     & ymin,ymax,       ! Min/Max igram values for the scan
     & infomat(mif,mns),! Information produced by real-time algorithm
     & Tdur,            ! DUR parameter as defined in original opus file
     & Tscan,           ! Time for 1 scan, calculated from NPT and VEL
     & Tturn            ! Turnaround time between fwd-rev scans

      character*40
     & DTCstr,          ! Character variable to hold DTC description
     & INSstr           ! Character variable to hold INS description

      data infovec/mif*0.0/
      idum = src_off  ! Avoid compiler warning (unused)
      idum = src_sun  ! Avoid compiler warning (unused)
      idum = src_mir  ! Avoid compiler warning (unused)
      idum = src_nir  ! Avoid compiler warning (unused)
      idum = sfm_script  ! Avoid compiler warning (unused)
      idum = sfm_idle ! Avoid compiler warning (unused)

      idum=i_aptval   ! Avoid compiler warnings (unused)
      idum=i_aqmcode  ! Avoid compiler warnings (unused)
      idum=i_bmscode  ! Avoid compiler warnings (unused)
      idum=i_dtccode  ! Avoid compiler warnings (unused)
      idum=i_inscode  ! Avoid compiler warnings (unused)

      idum=bms_caf2   ! Avoid compiler warnings (unused)
      idum=bms_kbr    ! Avoid compiler warnings (unused)
      idum=bms_quartz ! Avoid compiler warnings (unused)
      idum=bms_sica   ! Avoid compiler warnings (unused)

      idum=i_lpf      ! Avoid compiler warnings (unused)
      idum=i_mvd      ! Avoid compiler warnings (unused)
      idum=i_p2a      ! Avoid compiler warnings (unused)
      idum=i_p2k      ! Avoid compiler warnings (unused)
      idum=i_p2l      ! Avoid compiler warnings (unused)
      idum=i_p2r      ! Avoid compiler warnings (unused)
      idum=i_pgn      ! Avoid compiler warnings (unused)
      idum=i_pka      ! Avoid compiler warnings (unused)
      idum=i_pra      ! Avoid compiler warnings (unused)
      idum=i_prl      ! Avoid compiler warnings (unused)
      idum=i_pkl      ! Avoid compiler warnings (unused)
      idum=i_rsn      ! Avoid compiler warnings (unused)
      idum=i_sgna     ! Avoid compiler warnings (unused)
      idum=i_sgnb     ! Avoid compiler warnings (unused)
      idum=i_srccode  ! Avoid compiler warnings (unused)
      idum=i_sfmcode  ! Avoid compiler warnings (unused)
      idum=i_ssm      ! Avoid compiler warnings (unused)
      idum=i_ssp      ! Avoid compiler warnings (unused)
      idum=i_zff      ! Avoid compiler warnings (unused)
      idum=i_lwn      ! Avoid compiler warnings (unused)

      idum=i_dur      ! Avoid compiler warnings (unused)
      idum=i_foc      ! Avoid compiler warnings (unused)
      idum=i_hpf      ! Avoid compiler warnings (unused)

      idum=sfm_solar  ! Avoid compiler warnings (unused)
      idum=sfm_script ! Avoid compiler warnings (unused)
      idum=sfm_cell   ! Avoid compiler warnings (unused)
      idum=sfm_aeros  ! Avoid compiler warnings (unused)
      idum=aqm_sd     ! Avoid compiler warnings (unused)
      idum=aqm_sf     ! Avoid compiler warnings (unused)
      idum=aqm_sn     ! Avoid compiler warnings (unused)

c  Program execution starts here: display version string.
      write(*,'(a)')verstr

c  Initialize variables.
      errnum=0
      inpstat=0
      scancnt=0
      i4head(i_bfw)=0    ! We never save bad runs: BFW=BBW=0.
      i4head(i_bbw)=0 
      sivcflag=0
      dclevel=0.0d0
      fvsi_calc=0.0d0
      shbar=0.    ! Missing/impossible value
      sherr=0.    ! Missing/impossible value
      fpilha=0.   ! Missing/impossible value

c  Check that input file exists, if so open it.
      call getarg(1,infile)
      istat=lnbc(infile)
      if(istat.le.0) then
         stop 'Input filename missing from command line '
      endif

      inquire(file=infile,exist=filexist,iostat=inpstat)
      if(inpstat.ne.0) then
         errnum=-1
         write(*,'(2a)')'Error: inquire failed on input file ',infile
      elseif(filexist) then
         open(unit=lunr_in,file=infile,status='old',iostat=inpstat)
         if(inpstat.ne.0) then
            errnum=-1
            write(*,'(2a)')'Error: open failed on input file ',infile
         endif
      else
         errnum=-1
         write(*,'(2a)')'Error: cannot find input file ',infile
      endif

c   Check to see if the "Resistors" keyword is present.
      call read_input_line(lunr_in,errnum,inpstat,inpath)
      if((errnum.eq.0).and.(inpstat.ne.0)) then
         write(*,'(a)')'Error in input file: could not read igram path'
         errnum=-2
      elseif((errnum.eq.0).and.(index(inpath,"Resistors").eq.1)) then
         read(inpath,*,iostat=inpstat) inpstr,errnum
         if(inpstat.ne.0) then
            write(*,'(a)') 'Error in input file: bad resistor request'
            errnum=-2
         else
            write(*,*) 'Calling check_resistor'
c            call check_resistor_info(errnum)
            errnum=errnum+1000
         endif
      endif
      rewind(lunr_in)

c  Parse the top section of the input file containing general parameters.
      call parse_input_top(lunr_in,mdtc,errnum,inpstat,inpath,outpath,
     & igrmode,igrmpath,phmode,phpath,chan1,chan2,flimit,pattern,
     & srcindic,dtcigrm,dtcspec,outfmt,delimit,minmax,stlimstd,
     & stlimavg,ylimits,xcorlim,timecorr,mscan,fftlim,sivcfreq,
     & lsemode,pco_len,pco_thresh,proclim,verbose,run_start)

c  Initialize the memory of date to force reset of the run number.
      iyyywas=-1
      imwas=-1
      idwas=-1
c      memyear=-1
c      memmonth=-1
c      memday=-1
c      catslice=-1    ! Avoid catalog work-around on first entry.
      catslice=0
      nil=0
      nbad=0
c  Main program loop, over the scan sets listed in the input file.
c  Repeat when not at EOF, no errors.
      do while( (inpstat.eq.0) .and. (errnum.eq.0) )

c  Position input file to next slice entry or to end-of-file.
         call read_input_line(lunr_in,errnum,inpstat,inpstr)
         if(inpstat.ne.0) exit
         nil=nil+1

c  Reset "Number of Sample Scans" in case there is an error detected
c  before it can be read from the slice/igram headers.
         nss=0
         scanllt = 1
         scancnt = 0

c  Parse one line of the input file.
c         slicmem=catslice
         call substr(inpstr,cc,1,nsubstr)
         if(nsubstr.eq.5) then  ! slice_i2s
            read(inpstr,*,iostat=inpstat)
     &      catyear,catmonth,catday,catbatch,catslice
            inlen=lnbc(inpath)
            write(fpsfname,'(3i2.2,a1,i0,a1,a4,2a1)')
     &      mod(catyear,100),catmonth,catday,'.',catbatch,
     &      inpath(inlen:inlen),'scan',inpath(inlen:inlen),'b'
            if(inpstat.ne.0) then
               write(*,'(a)')'Format error in catalog input file'
               errnum=1
               inpstat=0
c            elseif(slicmem.eq.catslice) then
c               catslice=catslice+1
c               write(*,'(a)')
c     &         'ACHTUNG! catalog work-around activated: incremented slice'
            endif

c  See if we need to skip a run number: each line of the input file
c  starts with a FWD, and FWD runs must have an odd number.
            if(errnum.eq.0) then
               if(mod(runno,2).eq.0) runno=runno+1
            endif
         elseif(nsubstr.eq.16) then  ! opus-i2s
            lb=fnbc(inpstr)
            lf=fbc(inpstr(lb:))+lb
            fpsfname=inpstr(lb:lf-1)
            read(inpstr(lf:),*,iostat=inpstat)
     &      catyear,catmonth,catday,runno,
     &      infovec(7),infovec(8),infovec(9),    ! LAT, LON, ALT
     &      infovec(4),infovec(3),infovec(1),    ! TSC, PIM, HUM
     &      infovec(15),infovec(19),infovec(16), ! TOU, POU, HOU
     &      infovec(28), infovec(29)             ! SIA, fvsi (fractional)
            infovec(29)=infovec(29)*infovec(28)  ! convert fvsi to absolute(SIS)
            infovec(2)=infovec(4)                ! TLP
         elseif(nsubstr.eq.18) then  ! opus-i2s
            lb=fnbc(inpstr)
            lf=fbc(inpstr(lb:))+lb
            fpsfname=inpstr(lb:lf-1)
            read(inpstr(lf:),*,iostat=inpstat)
     &      catyear,catmonth,catday,runno,
     &      infovec(7),infovec(8),infovec(9),    ! LAT, LON, ALT
     &      infovec(4),infovec(3),infovec(1),    ! TSC, PIM, HUM
     &      infovec(15),infovec(19),infovec(16), ! TOU, POU, HOU
     &      infovec(28), infovec(29),            ! SIA, fvsi (fractional)
     &      infovec(10), infovec(13)             ! WSPD, WDIR  added DG 110220
            infovec(29)=infovec(29)*infovec(28)  ! convert fvsi to absolute(SIS)
            infovec(2)=infovec(4)                ! TLP
         else
            stop 'unrecognized info file format'
         endif
         iyyy=catyear
         im=catmonth
         id=catday

         write(*,*) 'fpsfname = ',fpsfname

         if(inpstat.ne.0) then
            write(*,'(a)')'Format error in catalog input file'
            errnum=1
            inpstat=0
         endif

         if(infovec(28).ne.0.0 .and. infovec(29).ne.0.0) then
            if(stlimavg.ne.0.0.and.stlimstd.ne.0) then
               if(infovec(28).lt.stlimavg.or.
     &         infovec(29).gt.stlimstd*infovec(28)) then
                  write(*,'(3a,f6.1,a,f6.2)') fpsfname(:lnbc(fpsfname)),
     &            ': Tracker intensity too low or stdev too ',
     &            'high: SIA = ', infovec(28),', SIS = ',infovec(29)
                  nbad=nbad+1
                  cycle
               endif
            endif
         endif

c  Obtain all IFS12X parameters by reading the header.
         if((errnum.eq.0).and.(proclim.ge.1)) then
            call get_igram_run_parameters(inpath,fpsfname,
     &      catslice,runno,verbose,mns,msl,mip,mch,mi4,mr8,
     &      chan1,chan2,errnum,nptvec,bpdata,nss,tpx,timvec,
     &      timsli,runsta,runend,i4head,r8head,DTCstr,INSstr)
         endif

c  Obtain auxiliary information produced by real-time software.
         if((errnum.eq.0).and.(proclim.ge.2)) then
            if(nsubstr.eq.5) then  ! slice-i2s
               call default_run_info(mif,mns,infomat)
               call get_run_info(inpath,catyear,catmonth,catday,
     &         catbatch,catslice,runno,verbose,mif,mns,nss,
     &        i4head(i_sfmcode),stlimavg,stlimstd,runsta,timvec,infomat)
            endif
         endif

c  Setup for FWD-only and FWD+REV runs.
         Tscan=0.0d0
         Tturn=0.0d0
         if(nsubstr.eq.5) then  ! slice-i2s
            tla_ext=3  ! tla_ext_full
            scanllt=nss
         else                 ! opus-i2s
            tla_ext=1  !tla_ext_min
            scanllt = 1
            scancnt = nss
            if((errnum.eq.0).and.(proclim.ge.1)) then
               if(i4head(i_aqmcode).eq.1 .and. nss.gt.1) then  ! aqm_sd
                  scanllt = 2
                  scancnt = nss/2
c     DG090401:  calculate scan time and turnaround time 
                  Tdur = r8head(i_dur)   ! Total duration time first-last point
                  Tscan = dble(tpx/2)/2.0D0/r8head(i_laserate)
     &                 /dble(i4head(i_ssm))*dble(i4head(i_ssp))
                  Tturn = (Tdur - nss*Tscan)/dble(nss-1)
                  r8head(i_dur) = Tdur - Tscan - Tturn
c     DG090401 end
                  if(mod(runno,2).ne.1) then
                     errnum=2
                     write(*,'(a)')'Error: FWD run # must be odd'
                  endif

                  if(mod(tpx,2).ne.0) then
                     errnum=2
                     write(*,'(3a,i0)')'Error:',
     &               fpsfname(1:lnbc(fpsfname)),' has odd TPX of ',tpx
                  else
                     tpx=tpx/2
                  endif

               endif   !  (i4head(i_aqmcode).eq.1 .and. nss.gt.1)
            endif   !  ((errnum.eq.0).and.(proclim.ge.1))
         endif   !  (nsubstr.ne.5)

c  Setup the OPUS header items that will remain unchanged throughout
c  this scan set.
         if((errnum.eq.0).and.(proclim.ge.1)) then
            if(nsubstr.eq.5) then ! slice-i2s
c  The running sample number (i_rsni, a "serial number") is set to
c  the starting slice number.
               i4head(i_rsn)=catslice
c  The scan duration (idur) is computed as the number of data points
c  divided by data rate.
               r8head(i_dur)=dble(tpx)/((2.d0*r8head(i_laserate)*
     &         dble(i4head(i_ssm)))/dble(i4head(i_ssp)))
            endif

c  FIXME: the zero-filling factor could be set to values greater than 2
c  if we allowed the FFT size to grow beyond the nearest power of two.
            i4head(i_zff)=2

c  Compute folding limits for this spectral domain.  We don't use the
c  information from the igram headers because it is not accurate enough.
c  We did make sure that LWN is correct, so we use that instead.
            minfold=0.d0
            maxfold=minfold+((r8head(i_lwn)*
     &      dble(i4head(i_ssm)))/dble(i4head(i_ssp)))
         endif

c  Now we process each scan.
c  Scans that were rejected or found in error have a negative 'runsta'.
         isi=1
         do iscan=1,scanllt        ! FWD then REV

c  Shift REV scan (iscan=2) into upper half of memory (for opus-i2s).
            if(nsubstr.ne.5) isi=1+(iscan-1)*mip/2

c  Apply site-specific time correction.
            timvec(iscan)=timvec(iscan)+(dble(timecorr)*3600.d0)


            if( nsubstr.gt.5 .or. runsta(iscan).gt.0 ) then  ! good opus-ipp scan

               write(*,*)' Good scan: iscan,nsubstr,runsta = ',
     &         iscan,nsubstr,runsta(iscan)

               if(nsubstr.eq.5) then ! slice-ipp

c  Copy the column of 'infomat' for this scan into 'infovec'.
                  do ia=1,mif
                     infovec(ia)=infomat(ia,iscan)
                  enddo

c  Set run direction info based on parity of run number (odd=fwd).
                  if(mod(runno,2).ne.0) then
                     i4head(i_gfw)=1
                  else
                     i4head(i_gfw)=0
                  endif
                  i4head(i_gbw)=1-i4head(i_gfw)

c  Now we process each channel from this one scan.
                  inlen=lnbc(inpath)
                  write(fpsfname,'(3i2.2,a1,i0,a1,a4,2a1)')
     &            mod(catyear,100),catmonth,catday,'.',catbatch,
     &            inpath(inlen:inlen),'scan',inpath(inlen:inlen),'b'

               endif  !   if(nsubstr.eq.5)

c   New section of code inserted Feb 26, 2010 by GCT.
c   Computes local scan date from OPUS header info (not catalog)
c   This local scan date is subsequently used by build_cit_name.
c   TIMVEC contains igram time info (UT) from OPUS header.
c   INFOVEC(8) is site longitude, needed to correct UT date to local.
               call julian(2000,1,1,j0)
               jd=j0+int((timvec(iscan)/86400)+infovec(8)/360)
               call caldat(jd,iyyy,im,id)
               if(iyyy.ne.iyyywas.or.im.ne.imwas.or.id.ne.idwas) then
                  write(*,*)'Detected date change within saveset'
                  if(nsubstr.eq.5) runno=run_start+mod(runno-1,2)
                  iyyywas=iyyy
                  imwas=im
                  idwas=id
               endif

               do ichan=chan2,chan1,-1  ! Si then InGaAs

c  Read interferogram(s)   (combined FWD+REV if AQM=SD and nss>1)
                  if((errnum.eq.0).and.(proclim.ge.3)) then
c                     write(*,*) ' Calling get_igram_data...'
                     call get_igram_data(inpath,fpsfname,catslice,
     &               runsta(iscan),runend(iscan),verbose,
     &               msl,mip,mch,nptvec,bpdata,ichan,nip,buf_igram)

c  Following 5 lines from JFB email of 2018-07-11 (counter --> nip)
                     if(nip.eq.(2*tpx) .and. nsubstr.eq.5) then
                        write(*,*)'SINGLE-SLICE MULTI-SCAN SLICE-I2S
     &  RUN - DOING FWD ONLY'
                        nip=tpx
                     endif

c  Find min/max igram values
                     ymin=buf_igram(1)
                     ymax=buf_igram(1)
                     do i=2,nip
                        if(buf_igram(i).lt.ymin) ymin=buf_igram(i)
                        if(buf_igram(i).gt.ymax) ymax=buf_igram(i)
                     end do

c  Reject scans if they contain igram values out of the allowed range:
c      Ymin <  MIN_thresh
c      Ymax >  MAX_thresh
c  ylimits(ichan,1) is the min allowed y-value (entry #17 of input file).
c  ylimits(ichan,2) is the max allowed y-value (entry #17 of input file).

c   Is Ymin below MIN_thresh?
                     if( ymin .lt. ylimits(ichan,1) ) then
                        write(*,*)'Scan rejected. Ymin < MIN_thresh',
     &                  catyear,catmonth,catday,runno+scancnt-1,iscan
                        write(*,*)'ichan,Ymin,MIN_thresh =',
     &                  ichan,ylimits(ichan,1),ylimits(ichan,2)
c                        runsta(iscan)=-101
                        errnum=+111
                     endif

c   Is Ymax above MAX_thresh?
                     if( ymax .gt. ylimits(ichan,2) ) then
                        write(*,*)'Scan rejected. Ymax > MAX_thresh',
     &                  catyear,catmonth,catday,runno+scancnt-1,iscan
                        write(*,*)'ichan,Ymax,MAX_thresh =',
     &                  ichan,ymax,ylimits(ichan,2)
c                        runsta(iscan)=-101
                        errnum=+112
                     endif

c  If file was recorded in SD mode and contains more than one scan,
c  divide NIP by two to separate FWD from REV runs. 
c  originally, FWD scans cover 1,nip and REV scans nip+1,2*nip
c FIXME: check that TPX and NIP match
                     if(nsubstr.ge.16) then  ! opus-ipp
                        if((errnum.eq.0).and.(scanllt.eq.2)) then
                           nip=nip/2
                           do ia=nip,1,-1
                              buf_igram(mip/2+ia)=buf_igram(nip+ia)
                           end do
                        endif
c  Now FWD scans still cover 1,nip  REV scans mip/2+1,mip/2+nip

c  Set run direction variable based on parity of run number (odd=fwd).
c  For reverse runs, get data back from temporary buf_igram.
c  DG090226 set TIM as start time for fwd and rev scans 

c  DG090522 The next line crashed Lauder runs with only 
c  FWD single scans and consecutive run numbers
c                 if(mod(runno,2).ne.0) then  !FWD scan
c  So replaced with:
                        if(scanllt.eq.1) then !single scans
                           if(i4head(i_gbw).eq.1) then ! SLICE-I2S REV scans
                              i4head(i_gfw)=0
c SLICE-I2S generated reverse scan: want buf_igram,
                              Tstart = timvec(iscan)+Tscan+Tturn
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
                              Tstart = timvec(iscan)
                           endif
                        else
                           if(mod(runno,2).ne.0) then  !FWD scan
                              i4head(i_gfw)=scancnt
                              Tstart = timvec(iscan)
                           else                        !REV scan
                              i4head(i_gfw)=0
                              do ia=nip,1,-1
                                 buf_igram(ia)=buf_igram(nip+ia)
                              end do
                              Tstart = timvec(iscan)+Tscan+Tturn
                           endif
                        endif
c                        i4head(i_gbw)=scancnt-i4head(i_gfw)
                        i4head(i_gbw)=1-i4head(i_gfw)   !  GCT kludge

                     else   ! slice-i2s
                        Tstart = timvec(iscan)   !  slice-i2s
                     endif   !   if(nsubstr.ge.16)

c  If enabled, save separated interferograms in the requested format.
                     if(igrmode.eq.1) then
                        call build_cit_name(igrmpath,pattern,dtcigrm,
     &                  srcindic,iyyy,im,id,ichan,i4head(i_srccode),
     &                  runno,filename)

                        call save_to_file(1,flimit,filename,ichan,
     &                  outfmt,verbose,progver,mip,mif,mi4,mr8,
     &                  buf_igram(isi),nip,nside,minfold,maxfold,
     &                  delimit,minmax,tpx,nshort,Tstart,i4head,r8head,
     &                  DTCstr,INSstr,sivcfreq(ichan),pco_len(ichan),
     &                  pco_thresh(ichan),izpd,sivcflag,dclevel,
     &                  fvsi_calc,zpa,frzpda,shbar,sherr,lsemode(ichan),
     &                  fpilha,infovec,tla_ext,fpsfname,errnum)
                     endif
                  endif

                  frsp=sngl((sivcfreq(ichan)-minfold)/(maxfold-minfold))
                  margin=nburst+pco_len(ichan)/2
c                  write(*,*) ' Calling siv_correction...'
                  call siv_correction(nip,margin,buf_igram(isi),
     &            smoo_igram,
     &            frsp,pinl,pinv,sivcflag,dclevel,fvsi_calc,frzpda)

c   Save SIV-corrected interferogram, if requested
                  if(igrmode.eq.2) then
                     call build_cit_name(igrmpath,pattern,dtcigrm,
     &               srcindic,iyyy,im,id,ichan,i4head(i_srccode),
     &               runno,filename)

                     call save_to_file(1,flimit,filename,ichan,outfmt,
     &               verbose,progver,mip,mif,mi4,mr8,buf_igram,nip,
     &               nside,minfold,maxfold,delimit,minmax,tpx,nshort,
     &               Tstart,i4head,r8head,DTCstr,INSstr,sivcfreq(ichan),
     &               pco_len(ichan),pco_thresh(ichan),izpd,sivcflag,
     &               dclevel,fvsi_calc,zpa,frzpda,shbar,sherr,
     &               lsemode(ichan),fpilha,infovec,tla_ext,
     &               fpsfname,errnum)
                  endif

c  Perform Fourier transform and save data in the requested format.
                  if((errnum.eq.0).and.(proclim.ge.4)) then
                     if(phmode.gt.0) then
                        call build_cit_name(phpath,pattern,dtcspec,
     &                  srcindic,iyyy,im,id,ichan,i4head(i_srccode),
     &                  runno,phasefull)
                        phasefull=phasefull(:lnbc(phasefull))//'.phs'
                     else
                        phasefull=''
                     endif

                  write(*,*)' calling i2s_processing...'

                     call i2s_processing(mip,fftlim(ichan),verbose,
     &               i4head(i_gfw),pinl,ichan,pco_len(ichan),
     &               pco_thresh(ichan),nburst,phasefull,
     &               buf_igram(isi),nip,nlong,nshort,nffthr,nside,
     &               izpd,zpa,shbar,sherr,lsemode(ichan),fpilha)
                  write(*,*)' called i2s_processing.'

                     call build_cit_name(outpath,pattern,dtcspec,
     &               srcindic,iyyy,im,id,ichan,i4head(i_srccode),
     &               runno,filename)

c  Save spectrum to file
                     call save_to_file(2,flimit,filename,ichan,outfmt,
     &               verbose,progver,mip,mif,mi4,mr8,buf_igram(isi),
     &               nffthr,nside,minfold,maxfold,delimit,minmax,
     &               nlong,nshort,Tstart,i4head,r8head,DTCstr,INSstr,
     &               sivcfreq(ichan),pco_len(ichan),pco_thresh(ichan),
     &               izpd,sivcflag,dclevel,fvsi_calc,zpa,frzpda,
     &               shbar,sherr,lsemode(ichan),fpilha,infovec,tla_ext,
     &               fpsfname,errnum)
                  endif     !  (errnum.eq.0).and.(proclim.ge.4)
               enddo      ! ichan=chan1,chan2
               if(errnum.eq.0) scancnt=scancnt+1
            else
               write(*,*)'Skipping bad scan: iscan,nsubstr,runsta = ',
     &         iscan,nsubstr,runsta(iscan)
            endif       ! (runsta(iscan).gt.0)

c  Increment the run number used for file names, unless an idle scan.
            if(runsta(iscan).ne.(-100)) runno=runno+1
         write(*,*)'--------- Loop over FWD/REV -------------'
         enddo          ! iscan=1,scanllt

c  Positive error codes indicate transient errors, specific to one set
c  of NSS scans.  The program can proceed to the next set of scans.
         if(errnum.gt.0) errnum=0

         write(*,*)'============ Loop over scan sets =============='
      enddo      ! while (inpstat.eq.0).and.(errnum.eq.0)

c  Close the input file.
      close(unit=lunr_in,iostat=inpstat)
      if(inpstat.ne.0) write(*,'(2a)')'Error: close failed on ',infile

      print *
      print '(2(I6,a))', nil, ' opus files in catalog, ',
     &                   nbad, ' excluded by SIA/SIS'

c  Terminate.
      stop
      end
