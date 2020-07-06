      subroutine save_opus(path,verbose,progver,mip,mif,mi4,mr8,
     & datype,ichan,buffer,counter,minfrq,maxfrq,minind,maxind,
     & nlong,nshort,time,i4head,r8head,DTCstr,INSstr,sivcfrq,
     & pco_leni,pco_threshi,
     & izpd,sivcflag,dclevel,fvsi_calc,zpa,frzpda,
     & shbar,sherr,lsemode,fpilha,snr,
     & infovec,tlalevel,fpsfname,errnum)
c
c  Input:
c    path      C*(*)  Path to OPUS output file
c    verbose   I*4    Level of verbosity for displayed messages
c    progver   R*8    Program version (date) to be stored in OPUS header
c    mip       I*4    Maximum number of input points
c    mif       I*4    Maximum number of info channels
c    mi4       I*4    Maximum number of I*4 items in file header
c    mr8       I*4    Maximum number of R*8 items in file header
c    datype    I*4    Data type (1=interferogram, 2=spectrum)
c    channel   I*4    Channel number (1=InGaAs=slave, 2=Si=master)
c    buffer(mip)R*4   Buffer containing interferogram/spectrum data
c    counter   I*4    Count of points in buffer
c    minfrq   R*8    Frequency of the first data point in buffer
c    maxfrq   R*8    Frequency one spectral point past the last one in buffer
c    minind    I*4    Point index of the lowest wavenumber to be saved
c    maxind    I*4    Point index of the highest wavenumber to be saved
c    nlong     I*4    Number of points in one FWD or REV scan
c    nshort    I*4    Number of points in one FWD or REV scan
c    time      R*8    Time in seconds since 1-Jan-2000
c    i4head(mi4)I*4   Vector holding the I*4 header items
c    r8head(mr8)R*8   Vector holding the R*8 header items
c    izpd      I*4    Point index of zero path difference
c    infovec(mif)R*8  Information produced by real-time algorithm
c    tlalevel  I*4    Level of non-Bruker header items (three-letter acronyms)
c
c  Input/Output:
c    errnum    I*4    Program error status
c
      implicit none

      integer*4
     & idum,lp,
     & verbose,    ! Subroutine input argument (see above)
     & mip,        ! Subroutine input argument (see above)
     & mif,        ! Subroutine input argument (see above)
     & mi4,        ! Subroutine input argument (see above)
     & mr8,        ! Subroutine input argument (see above)
     & datype,     ! Subroutine input argument (see above)
     & ichan,    ! Subroutine input argument (see above)
     & counter,    ! Subroutine input argument (see above)
     & minind,     ! Subroutine input argument (see above)
     & maxind,     ! Subroutine input argument (see above)
     & nlong,      ! Subroutine input argument (see above)
     & nshort,     ! Subroutine input argument (see above)
     & i4head(mi4),! Subroutine input argument (see above)
     & sivcflag,    !
     & lsemode,    ! Laser sampling error type
     & pco_leni,
     & izpd,       ! Subroutine input argument (see above)
     & tlalevel,   ! Subroutine input argument (see above)
     & errnum,     ! Subroutine input/output argument (see above)
     & iend,       ! Endianness of computer
     & bigendian,  ! Named constant for big endian detection
     & lnbc,       ! Integer function Last Non-Blank Character in string
     & indexa,     ! General loop index
     & luna,       ! Logical Unit Number for output file
     & ierr,       ! Error status returned by file I/O
     & ndb,        ! Current number of directory blocks
     & pointr,     ! Main byte pointer into the file header
     & bytstart,   ! Starting byte pointer for each parameter block
     & length,     ! Length in I*4 of each parameter block
     & headir,     ! Location in header block of pointer to directory
     & dirpnt,     ! Byte offset for current directory entry
     & ndbpnt,     ! Location in header block of size of directory
     & idtagd,     ! Block type for the data
     & idtags,     ! Block type for the data status parameters
     & year,       ! Date of the first data point in this scan
     & month,      ! Date of the first data point in this scan
     & day,        ! Date of the first data point in this scan
     & hour,       ! Time of the first data point in this scan
     & minute,     ! Time of the first data point in this scan
     & second,     ! Time of the first data point in this scan
     & millisec,   ! Time of the first data point in this scan
     & imxy,       ! Index of maximum Y value
     & imny        ! Index of minimum Y value

      real*4
     & buffer(mip), ! Subroutine input argument (see above)
     & snr,
     & shbar,       ! LSE
     & sherr,       ! LSU
     & fpilha       ! LSF

      real*8
     & progver,     ! Subroutine input argument (see above)
     & minfrq,      ! Subroutine input argument (see above)
     & maxfrq,      ! Subroutine input argument (see above)
     & time,        ! Subroutine input argument (see above)
     & r8head(mr8), ! Subroutine input argument (see above)
     & infovec(mif),!Subroutine input argument (see above)
     & res,         ! Spectral resolution, in cm-1
     & phr,         ! Phase resolution, in cm-1
     & sivcfrq,     ! SIV-Correction Frequency (cm-1)
     & pco_threshi,
     & dclevel,
     & fvsi_calc,
     & zpa,         ! ZPD interferogram amplitude (phase-corrected)
     & frzpda,
     & fxv,         ! First X value, in point index (igrams) or cm-1 (spectra)
     & lxv,         ! Last X value, in point index (igrams) or cm-1 (spectra)
     & mxy,         ! Maximum Y value
     & mny,         ! Minimum Y value
     & r8val,       ! Temporary variable used during min/max calculations
     & hfl          ! High folding limit, in cm-1

      character
     & fpsfname*(*),
     & path*(*),    ! Subroutine input argument (see above)
     & stringa*128, ! String used to format header items
     & DTCstr*(*),  ! Detector description
     & INSstr*(*)   ! Detector description

      logical*4
c     & filexist,   ! Keeps track of file existence
     & swapped      ! Indicates if we byte-swapped the data vector

      integer*4
c     & tla_ext_none, ! Use only the Bruker-defined header items
     & tla_ext_min,  ! Add those that apply to all TCCON sites
     & tla_ext_pgr,  ! Also add the Pre-amp Gain Resistor value
     & tla_ext_full  ! Complete set from automated system (Zeno, ST, dome,...)
c      parameter (tla_ext_none = 0)
      parameter (tla_ext_min  = 1)
      parameter (tla_ext_pgr  = 2)
      parameter (tla_ext_full = 3)

      parameter (bigendian=1)
      parameter (luna=19)

      include 'opus_constants.inc'
      include 'header_indices.inc'

      integer*1 header(mhl)! Work space to build the OPUS file header

      idum=i_ssm        !  Avoid compiler warnings (unused parameter)
      idum=i_ssp        !  Avoid compiler warnings (unused parameter)
      idum=i_srccode    !  Avoid compiler warnings (unused parameter)
      idum=i_sfmcode    !  Avoid compiler warnings (unused parameter)
      idum=i_dtccode    !  Avoid compiler warnings (unused parameter)
      idum=i_pgn        !  Avoid compiler warnings (unused parameter)
      idum=i_sgna       !  Avoid compiler warnings (unused parameter)
      idum=i_sgnb       !  Avoid compiler warnings (unused parameter)
      idum=i_hpf        !  Avoid compiler warnings (unused parameter)
      idum=i_gfw        !  Avoid compiler warnings (unused parameter)
      idum=i_gbw        !  Avoid compiler warnings (unused parameter)
      idum=i_bfw        !  Avoid compiler warnings (unused parameter)
      idum=i_bbw        !  Avoid compiler warnings (unused parameter)
      idum=i_rsn        !  Avoid compiler warnings (unused parameter)
      idum=i_pka        !  Avoid compiler warnings (unused parameter)
      idum=i_pkl        !  Avoid compiler warnings (unused parameter)
      idum=i_pra        !  Avoid compiler warnings (unused parameter)
      idum=i_prl        !  Avoid compiler warnings (unused parameter)
      idum=i_p2a        !  Avoid compiler warnings (unused parameter)
      idum=i_p2l        !  Avoid compiler warnings (unused parameter)
      idum=i_p2r        !  Avoid compiler warnings (unused parameter)
      idum=i_p2k        !  Avoid compiler warnings (unused parameter)
      idum=i_zff        !  Avoid compiler warnings (unused parameter)
      idum=i_inscode    !  Avoid compiler warnings (unused parameter)
      idum=i_bmscode    !  Avoid compiler warnings (unused parameter)
      idum=i_aqmcode    !  Avoid compiler warnings (unused parameter)
      idum=i_laserate   !  Avoid compiler warnings (unused parameter)
      idum=i_lwn        !  Avoid compiler warnings (unused parameter)
      idum=i_foc        !  Avoid compiler warnings (unused parameter)
      idum=i_aptval     !  Avoid compiler warnings (unused parameter)
      idum=i_lpf        !  Avoid compiler warnings (unused parameter)
      idum=i_dur        !  Avoid compiler warnings (unused parameter)
      idum=i_mvd        !  Avoid compiler warnings (unused parameter)
c GCT   idum=tla_ext_none !  Avoid compiler warnings (unused parameter)
c GCT   idum=tla_ext_min  !  Avoid compiler warnings (unused parameter)
c GCT   idum=tla_ext_pgr  !  Avoid compiler warnings (unused parameter)
c GCT   idum=tla_ext_full !  Avoid compiler warnings (unused parameter)
      idum=src_sun      !  Avoid compiler warnings (unused parameter)
      idum=src_nir      !  Avoid compiler warnings (unused parameter)
      idum=src_mir      !  Avoid compiler warnings (unused parameter)
      idum=src_off      !  Avoid compiler warnings (unused parameter)
      idum=sfm_solar    !  Avoid compiler warnings (unused parameter)
      idum=sfm_cell     !  Avoid compiler warnings (unused parameter)
      idum=sfm_idle     !  Avoid compiler warnings (unused parameter)
      idum=sfm_script   !  Avoid compiler warnings (unused parameter)
      idum=sfm_aeros    !  Avoid compiler warnings (unused parameter)
      idum=bms_caf2     !  Avoid compiler warnings (unused parameter)
      idum=bms_kbr      !  Avoid compiler warnings (unused parameter)
      idum=bms_quartz   !  Avoid compiler warnings (unused parameter)
      idum=bms_sica     !  Avoid compiler warnings (unused parameter)
      idum=aqm_sd       !  Avoid compiler warnings (unused parameter)
      idum=aqm_sf       !  Avoid compiler warnings (unused parameter)
      idum=aqm_sn       !  Avoid compiler warnings (unused parameter)
      idum=i_inscode    ! Avoid compiler warnings (unused parameter)  
      idum=typ_i4       ! Avoid compiler warnings (unused parameter)
      idum=typ_r8       ! Avoid compiler warnings (unused parameter)
c
c  Initialize variables.
c
      swapped=.false. 
c
c  The first part of this subroutine computes the values of some of the
c  header parameters.
c
c  Initialize identification tag of the data block.
c
c      write(*,*)' Inside save_opus.f'
      if(errnum.eq.0) then
         idtagd=dbbampl+dbbsamp      ! File has CPLX=amplitude and TYP=sample
         if(datype.eq.1) then
            idtagd=idtagd+dbbigrm     ! Either DATA=interferogram
         elseif(datype.eq.2) then
            idtagd=idtagd+dbbspec     ! or DATA=spectrum, undefined Y-units
         else
            errnum=-10
            write(*,'(a)')'Error in save_opus: invalid data type'
         endif
         if(ichan.eq.1) then
            idtagd=idtagd+dbbslav     ! Slave detector channel bit
         elseif(ichan.ne.2) then
            errnum=-10
            write(*,'(a)')'Error in save_opus: invalid channel'
         endif
      endif   ! (errnum.eq.0)
c
c  Initialize identification tag of the data status parameter block.
c
      if(errnum.eq.0) then
         idtags=idtagd+dbbdstat      ! PARM=data status parameter
c
c  Select the file limits based on data type.  For interferograms,
c  the X-scale is in units of data points, whereas for spectra the
c  units are wavenumbers.
c
         if(datype.eq.1) then
            fxv=dble(minind-1)
            lxv=dble(maxind-1)
         else
            fxv=minfrq+((dble(minind-1)*(maxfrq-minfrq))/dble(counter))
            lxv=minfrq+((dble(maxind-1)*(maxfrq-minfrq))/dble(counter))
         endif
c
c  Compute min/max of data buffer, needed for OPUS header.  Also, the
c  ZPD position is not yet known when saving interferograms, so set it
c  to the position of the extremum in that case.
c
         mxy=dble(buffer(minind))
         imxy=minind
         mny=mxy
         imny=imxy
         do indexa=minind+1,maxind
            r8val=dble(buffer(indexa))
            if(r8val.gt.mxy) then
               mxy=r8val
               imxy=indexa
            elseif(r8val.lt.mny) then
               mny=r8val
               imny=indexa
            endif
         enddo
         if(datype.eq.1) then  ! igram
            if(dabs(mxy).gt.dabs(mny)) then
               izpd=imxy
            else
               izpd=imny
            endif
            if(i4head(i_gbw).gt.i4head(i_gfw)) then   ! Reverse run
               izpd=nlong
            endif 
         endif
c
c  Compute the high folding limit.
c
         hfl=r8head(i_lwn)*dble(i4head(i_ssm))/dble(i4head(i_ssp))
c
c  Compute spectrum and phase resolution, using the stupid Bruker definition
c  of 0.9/L.  By the Nyquist theorem, 1/(2*hfl) is the point spacing, hence
c  L = npts/(2*hfl).  And finally, res = 0.9/L = 0.9*2*hfl/npts
c
c         res=1.8d0*hfl/dble(tpx-izpd)
         
         res=1.8d0*hfl/dble(nlong-pco_leni/2)  !  GCT 2014-04-15
         phr=1.8d0*hfl/dble(pco_leni/2)
c
c  Convert from seconds since 1-Jan-2000 (stardate) to date and time.
c
         call stardate_to_date (time,
     &   year,month,day,hour,minute,second,millisec)

         if(verbose.ge.4) then
            write(*,'(i4,a1,5(i2.2,a1),i3.3)') 
     &      year,'-',month,'-',day,' ',hour,
     &      ':',minute,':',second,'.',millisec
         endif
c
      endif   ! (errnum.eq.0)
c
c  Start building the OPUS file structure.  As much as possible, the
c  structure is built from the start of the file.  When we get to an
c  item whose value depends on subsequent items, we keep a pointer to
c  it and move on.  But first make sure there were no errors so far.
c
      if(errnum.eq.0) then
c
c  Initialize variable tracking computer endianness.
c
         call getendian(iend)
c
c  Initialize the main pointer into the byte array used to build the
c  header+directory+parameter blocks.  Also initialize the directory
c  block counter.
c
         pointr=0
         ndb=0
c
c  The header block is first.  It contains five items: magic number,
c  file format version, pointer to directory, maximum size of directory,
c  and current size of directory.  We don't know the value of the
c  'pointer to directory' and 'current size of directory' at this time.
c
         call put_i4_byte(magic,iend,mhl,errnum,pointr,header)
         call put_r8_byte(filevers,iend,mhl,errnum,pointr,header)
         headir=pointr
         call put_i4_byte(0,iend,mhl,errnum,pointr,header)
         call put_i4_byte(mdb,iend,mhl,errnum,pointr,header)
         ndbpnt=pointr
         call put_i4_byte(0,iend,mhl,errnum,pointr,header)
c
c  Now that we've reached the end of the header block, we know the
c  value of the pointer to the directory and we can put it back
c  into the header block.
c
         dirpnt=pointr
         call put_i4_byte(dirpnt,iend,mhl,errnum,headir,header)
c
c  Start building the directory structure.  The first entry in the
c  directory points to itself.  Don't ask.  The block type 'dbbdir'
c  has value 13*(2^10) = 13312.
c
         call put_i4_byte(dbbdir,iend,mhl,errnum,pointr,header)
         call put_i4_byte(3*mdb,iend,mhl,errnum,pointr,header)
         call put_i4_byte(dirpnt,iend,mhl,errnum,pointr,header)
         ndb=ndb+1
c
c  We will build the rest of the directory as we go: save the current
c  unfilled location into 'dirpnt'.  Then zero out the rest of the
c  directory which advances 'pointr' to the first parameter block.
c
         dirpnt=pointr
         do indexa=1,3*(mdb-ndb)
            call put_i4_byte(0,iend,mhl,errnum,pointr,header)
         enddo
c
c  Create instrument status parameter block.  The block type 'dbbinstr'
c  has value 2*(2^4) = 32.
c
         bytstart=pointr
         call pk_r8_opus('DUR',r8head(i_dur),iend,mhl,errnum,
     &   pointr,header)
         call pk_r8_opus('MVD',r8head(i_mvd),iend,mhl,errnum,
     &   pointr,header)

         if(ichan.eq.1) then
            call pk_i4_opus('PKA',i4head(i_p2a),iend,mhl,errnum,
     &      pointr,header)
            call pk_i4_opus('PKL',i4head(i_p2l),iend,mhl,errnum,
     &      pointr,header)
            call pk_i4_opus('PRA',i4head(i_p2r),iend,mhl,errnum,
     &      pointr,header)
            call pk_i4_opus('PRL',i4head(i_p2k),iend,mhl,errnum,
     &      pointr,header)
         else
            call pk_i4_opus('PKA',i4head(i_pka),iend,mhl,errnum,
     &      pointr,header)
            call pk_i4_opus('PKL',i4head(i_pkl),iend,mhl,errnum,
     &      pointr,header)
            call pk_i4_opus('PRA',i4head(i_pra),iend,mhl,errnum,
     &      pointr,header)
            call pk_i4_opus('PRL',i4head(i_prl),iend,mhl,errnum,
     &      pointr,header)
         endif
         call pk_r8_opus('HFL',hfl,iend,mhl,errnum,pointr,header)
         call pk_r8_opus('LFL',0.d0,iend,mhl,errnum,pointr,header)
         call pk_r8_opus('LWN',r8head(i_lwn),iend,mhl,errnum,
     &   pointr,header)
         call pk_i4_opus('SSP',i4head(i_ssp),iend,mhl,errnum,
     &   pointr,header)
         call pk_i4_opus('SSM',i4head(i_ssm),iend,mhl,errnum,
     &   pointr,header)

c     DG090401: start
c     Copy instrument description direct from input igm file
         call pk_st_opus('INS',INSstr,typ_string,iend,
     &   mhl,errnum,pointr,header)

c         if(i4head(i_inscode).eq.ins_120) then
c            call pk_st_opus('INS','IFS120HR',typ_string,iend,
c     &       mhl,errnum,pointr,header)
c         elseif(i4head(i_inscode).eq.ins_120m) then
c            call pk_st_opus('INS','IFS120M',typ_string,iend,
c     &      mhl,errnum,pointr,header)
c         elseif(i4head(i_inscode).eq.ins_125) then
c            call pk_st_opus('INS','IFS125HR',typ_string,iend,
c     &      mhl,errnum,pointr,header)
c         elseif(i4head(i_inscode).eq.ins_66s) then
c            call pk_st_opus('INS','IFS66/S',typ_string,iend,
c     &      mhl,errnum,pointr,header)
c         endif

         call pk_r8_opus('FOC',r8head(i_foc),iend,mhl,errnum,
     &   pointr,header)
         call pk_i4_opus('ASS',i4head(i_gfw)+i4head(i_gbw),iend,mhl,
     &   errnum,pointr,header)
         call pk_i4_opus('GFW',i4head(i_gfw),iend,mhl,errnum,
     &   pointr,header)
         call pk_i4_opus('GBW',i4head(i_gbw),iend,mhl,errnum,
     &   pointr,header)
         call pk_i4_opus('BFW',i4head(i_bfw),iend,mhl,errnum,
     &   pointr,header)
         call pk_i4_opus('BBW',i4head(i_bbw),iend,mhl,errnum,
     &   pointr,header)
         call pk_st_opus('IGN',fpsfname,typ_enum,iend,mhl,errnum,
     &   pointr,header)
         call pk_i4_opus('RSN',i4head(i_rsn),iend,mhl,errnum,
     &   pointr,header)
         call pk_r8_opus('HFF',hfl,iend,mhl,errnum,pointr,header)
         call pk_r8_opus('LFF',0.d0,iend,mhl,errnum,pointr,header)
         call pk_r8_opus('HUM',infovec(1),iend,mhl,errnum,
     &   pointr,header)
         call pk_r8_opus('TLP',infovec(2),iend,mhl,errnum,
     &   pointr,header)
         call pk_r8_opus('PIM',infovec(3),iend,mhl,errnum,
     &   pointr,header)
         call pk_r8_opus('TSC',infovec(4),iend,mhl,errnum,
     &   pointr,header)
         if(tlalevel.ge.tla_ext_full) then
            call pk_r8_opus('IDA',infovec(33),iend,mhl,errnum,
     &      pointr,header)

            if(nint(infovec(35)).eq.0) then
               call pk_st_opus('ISS','OK',typ_string,iend,mhl,
     &         errnum,pointr,header)
            elseif(nint(infovec(35)).eq.1) then
               call pk_st_opus('ISS','Stress',typ_string,iend,mhl,
     &         errnum,pointr,header)
            elseif(nint(infovec(35)).eq.2) then
               call pk_st_opus('ISS','Error',typ_string,iend,mhl,
     &         errnum,pointr,header)
            elseif(nint(infovec(35)).eq.3) then
               call pk_st_opus('ISS','Reset',typ_string,iend,mhl,
     &         errnum,pointr,header)
            else
               call pk_st_opus('ISS','Unknown',typ_string,iend,mhl,
     &         errnum,pointr,header)
            endif
         endif

         call end_opus_prm(mhl,errnum,pointr,header)
         length=(pointr-bytstart)/4
         call put_i4_byte(dbbinstr,iend,mhl,errnum,dirpnt,header)
         call put_i4_byte(length,iend,mhl,errnum,dirpnt,header)
         call put_i4_byte(bytstart,iend,mhl,errnum,dirpnt,header)
         ndb=ndb+1
c
c  Create standard acquisition parameter block.  The block type 'dbbaqpar'
c  has value 3*(2^4) = 48.
c
         bytstart=pointr

         if(i4head(i_aqmcode).eq.aqm_sd) then
            call pk_st_opus('AQM','SD',typ_enum,iend,mhl,errnum,
     &      pointr,header)
         elseif(i4head(i_aqmcode).eq.aqm_sf) then
            call pk_st_opus('AQM','SF',typ_enum,iend,mhl,errnum,
     &      pointr,header)
         elseif(i4head(i_aqmcode).eq.aqm_sn) then
            call pk_st_opus('AQM','SN',typ_enum,iend,mhl,errnum,
     &      pointr,header)
         endif

         if(ichan.eq.1) then
            indexa=i4head(i_sgna)
         else
            indexa=i4head(i_sgnb)
         endif
         if((indexa.ge.0).and.(indexa.le.9)) then
            write(stringa,'(i1)')indexa
         else
            write(stringa,'(i2)')indexa
         endif
         call pk_st_opus('SGN',stringa,typ_enum,iend,mhl,errnum,
     &   pointr,header)

         call pk_r8_opus('RES',res,iend,mhl,errnum,pointr,header)
         if(tlalevel.ge.tla_ext_min) then
c            write(*,*)'save_opus: OPS: ',nshort
c            write(*,*)'save_opus: OPL: ',nlong
            call pk_r8_opus('OPL',dble(nlong)/(2.d0*hfl),iend,mhl,
     &      errnum,pointr,header)
            call pk_r8_opus('OPS',dble(nshort)/(2.d0*hfl),iend,mhl,
     &      errnum,pointr,header)
         endif
         call pk_i4_opus('NSS',
     &   i4head(i_gfw)+i4head(i_gbw)+i4head(i_bfw)+i4head(i_bbw),
     &   iend,mhl,errnum,pointr,header)
         call end_opus_prm(mhl,errnum,pointr,header)
         length=(pointr-bytstart)/4
         call put_i4_byte(dbbaqpar,iend,mhl,errnum,dirpnt,header)
         call put_i4_byte(length,iend,mhl,errnum,dirpnt,header)
         call put_i4_byte(bytstart,iend,mhl,errnum,dirpnt,header)
         ndb=ndb+1
c
c  Create Fourier Transform (FT) parameter block if we are saving a
c  spectrum.  The block type 'dbbftpar' has value 4*(2^4) = 64.
c
         if(datype.eq.2) then
            bytstart=pointr
            call pk_st_opus('APF','BX',typ_enum,iend,mhl,errnum,
     &      pointr,header)
            call pk_r8_opus('PHR',phr,iend,mhl,errnum,pointr,
     &      header)
            call pk_r8_opus('SVF',sivcfrq,iend,mhl,errnum,pointr,
     &      header)
            call pk_i4_opus('PCL',pco_leni,iend,mhl,errnum,pointr,
     &      header)
            call pk_r8_opus('PCT',pco_threshi,iend,mhl,errnum,pointr,
     &      header)
            call pk_i4_opus('ZFF',i4head(i_zff),iend,mhl,errnum,
     &      pointr,header)

            if(tlalevel.ge.tla_ext_min) then
c In the OPUS format, the point indices are 0-based, hence "izpd-1" 
               call pk_i4_opus('ZPL',izpd-1,iend,mhl,errnum,pointr,
     &         header)
               call pk_i4_opus('SCF',sivcflag,iend,mhl,errnum,pointr,
     &         header)
               call pk_r8_opus('DCL',dclevel,iend,mhl,errnum,pointr,
     &         header)
               call pk_r8_opus('ZPA',zpa,iend,mhl,errnum,pointr,
     &         header)
               call pk_r8_opus('DIP',frzpda,iend,mhl,errnum,pointr,
     &         header)
               call pk_r8_opus('VER',progver,iend,mhl,errnum,pointr,
     &         header)
               call pk_r8_opus('VDC',fvsi_calc,iend,mhl,errnum,pointr,
     &         header)
               call pk_i4_opus('LST',lsemode,iend,mhl,errnum,pointr,
     &         header)
               call pk_r8_opus('LSE',dble(shbar),iend,mhl,errnum,pointr,
     &         header)
               call pk_r8_opus('LSU',dble(sherr),iend,mhl,errnum,pointr,
     &         header)
               call pk_r8_opus('LSF',dble(fpilha),iend,mhl,errnum,
     &         pointr,header)
               call pk_r8_opus('SNR',dble(snr),iend,mhl,errnum,
     &         pointr,header)
            endif
            call end_opus_prm(mhl,errnum,pointr,header)
            length=(pointr-bytstart)/4
            call put_i4_byte(dbbftpar,iend,mhl,errnum,dirpnt,header)
            call put_i4_byte(length,iend,mhl,errnum,dirpnt,header)
            call put_i4_byte(bytstart,iend,mhl,errnum,dirpnt,header)
            ndb=ndb+1
         endif
c
c  Create optics parameter block.  The block type 'dbboptpar' has value
c  6*(2^4) = 96.
c
         bytstart=pointr
         if(i4head(i_srccode).eq.src_sun) then
            call pk_st_opus('SRC','Sun Tracker',typ_senum,iend,
     &      mhl,errnum,pointr,header)
         elseif(i4head(i_srccode).eq.src_nir) then
            call pk_st_opus('SRC','NIR',typ_senum,iend,
     &      mhl,errnum,pointr,header)
         elseif(i4head(i_srccode).eq.src_mir) then
            call pk_st_opus('SRC','MIR',typ_senum,iend,
     &      mhl,errnum,pointr,header)
         elseif(i4head(i_srccode).eq.src_off) then
            call pk_st_opus('SRC','Off All',typ_senum,iend,
     &      mhl,errnum,pointr,header)
         endif

         if(r8head(i_aptval).le.(9.95d0)) then
            write(stringa,'(f4.2,a3)')r8head(i_aptval),' mm'
         else
            write(stringa,'(f5.2,a3)')r8head(i_aptval),' mm'
         endif
         call pk_st_opus('APT',stringa,typ_senum,
     &   iend,mhl,errnum,pointr,header)
         call pk_r8_opus('FOV',1000.d0*r8head(i_aptval)/r8head(i_foc),
     &   iend,mhl,errnum,pointr,header)

         if(i4head(i_bmscode).eq.bms_caf2) then
            call pk_st_opus('BMS','CaF2',typ_senum,iend,
     &      mhl,errnum,pointr,header)
         elseif(i4head(i_bmscode).eq.bms_quartz) then
            call pk_st_opus('BMS','Quartz',typ_senum,iend,
     &      mhl,errnum,pointr,header)
         elseif(i4head(i_bmscode).eq.bms_sica) then
            call pk_st_opus('BMS','SiCa',typ_senum,iend,
     &      mhl,errnum,pointr,header)
         endif

         if(r8head(i_laserate).le.(9950.d0)) then
            write(stringa,'(f3.1,a4)')r8head(i_laserate)/1000.d0,' kHz'
         else
            write(stringa,'(f4.1,a4)')r8head(i_laserate)/1000.d0,' kHz'
         endif
         call pk_st_opus('VEL',stringa,typ_enum,iend,
     &   mhl,errnum,pointr,header)

c FIXME: separate data channel from detector and use these:
c        parameter (dtc_ingaas = 1)
c        parameter (dtc_si     = 2)
c
c     DG090401: start
c     Copy detector description direct from input igm file
         call pk_st_opus('DTC',DTCstr,typ_senum,iend,
     &   mhl,errnum,pointr,header)

c         if((i4head(i_dtccode).eq.1).and.(ichan.eq.1)) then
c            call pk_st_opus('DTC','InGaAs-AC',typ_senum,iend,
c     &      mhl,errnum,pointr,header)
c         elseif((i4head(i_dtccode).eq.1).and.(ichan.eq.2)) then
c            call pk_st_opus('DTC','Si-AC',typ_senum,iend,
c     &      mhl,errnum,pointr,header)
c         elseif((i4head(i_dtccode).eq.2).and.(ichan.eq.1)) then
c            call pk_st_opus('DTC','InGaAs-DC',typ_senum,iend,
c     &      mhl,errnum,pointr,header)
c         elseif((i4head(i_dtccode).eq.2).and.(ichan.eq.2)) then
c            call pk_st_opus('DTC','Si-DC',typ_senum,iend,
c     &      mhl,errnum,pointr,header)
c         elseif(i4head(i_dtccode).eq.dtc_insb) then
c            call pk_st_opus('DTC','InSb-AC',typ_senum,iend,
c     &      mhl,errnum,pointr,header)
c         elseif(i4head(i_dtccode).eq.dtc_bolom) then
c            call pk_st_opus('DTC','Bolometer',typ_senum,iend,
c     &      mhl,errnum,pointr,header)
cc FIXME         elseif(i4head(i_dtccode).eq.4) then
cc FIXME          call pk_st_opus('DTC','InSb-DC',typ_senum,iend,
cc FIXME           mhl,errnum,pointr,header)
c         endif

         if(tlalevel.ge.tla_ext_pgr) then
            if(ichan.eq.1) then                       ! InGaAs
               call pk_i4_opus('PGR',nint(infovec(5)),iend,mhl,errnum,
     &         pointr,header)
            elseif(i4head(i_dtccode).le.2) then          ! Si
               call pk_i4_opus('PGR',nint(infovec(6)),iend,mhl,errnum,
     &         pointr,header)
            elseif(i4head(i_pgn).eq.0) then             ! InSb, gain A
               call pk_i4_opus('PGR',1500,iend,mhl,errnum,pointr,header)
            elseif(i4head(i_pgn).eq.1) then             ! InSb, gain B
               call pk_i4_opus('PGR',3600,iend,mhl,errnum,pointr,header)
            elseif(i4head(i_pgn).eq.2) then             ! InSb, gain C
               call pk_i4_opus('PGR',25000,iend,mhl,errnum,pointr,
     &         header)
            elseif(i4head(i_pgn).eq.3) then             ! InSb, gain D = ref
               call pk_i4_opus('PGR',10000,iend,mhl,errnum,pointr,
     &         header)
            endif
         endif

         if(i4head(i_hpf).ge.0) then
            write(stringa,'(i1)')i4head(i_hpf)
            call pk_st_opus('HPF',stringa,typ_enum,iend,
     &      mhl,errnum,pointr,header)
         endif

         if(r8head(i_lpf).le.(9.95d0)) then
            write(stringa,'(f3.1)')r8head(i_lpf)
         else
            write(stringa,'(f4.1)')r8head(i_lpf)
         endif
         call pk_st_opus('LPF',stringa,typ_enum,iend,
     &    mhl,errnum,pointr,header)

         call end_opus_prm(mhl,errnum,pointr,header)
         length=(pointr-bytstart)/4
         call put_i4_byte(dbboptpar,iend,mhl,errnum,dirpnt,header)
         call put_i4_byte(length,iend,mhl,errnum,dirpnt,header)
         call put_i4_byte(bytstart,iend,mhl,errnum,dirpnt,header)
         ndb=ndb+1
c
c  Create sample origin parameter block.  The block type 'dbborgpar' has value
c  10*(2^4) = 160.
c
         bytstart=pointr
         if(i4head(i_sfmcode).eq.sfm_solar) then
            call pk_st_opus('SFM','Solar',typ_string,iend,
     &      mhl,errnum,pointr,header)
         elseif(i4head(i_sfmcode).eq.sfm_cell) then
            call pk_st_opus('SFM','Cell',typ_string,iend,
     &      mhl,errnum,pointr,header)
         elseif(i4head(i_sfmcode).eq.sfm_idle) then
            call pk_st_opus('SFM','Short Scan for inactivity',
     &      typ_string,iend,mhl,errnum,pointr,header)
         elseif(i4head(i_sfmcode).eq.sfm_script) then
            call pk_st_opus('SFM','scripting',typ_string,iend,
     &      mhl,errnum,pointr,header)
         elseif(i4head(i_sfmcode).eq.sfm_aeros) then
            call pk_st_opus('SFM','Aerosol',typ_string,iend,
     &      mhl,errnum,pointr,header)
         endif

         if(tlalevel.ge.tla_ext_min) then
            call pk_r8_opus('LAT',infovec(7),iend,mhl,errnum,
     &      pointr,header)
            call pk_r8_opus('LON',infovec(8),iend,mhl,errnum,
     &      pointr,header)
            call pk_r8_opus('ALT',infovec(9),iend,mhl,errnum,
     &      pointr,header)
            call pk_r8_opus('TOU',infovec(15),iend,mhl,errnum,
     &      pointr,header)
            call pk_r8_opus('POU',infovec(19),iend,mhl,errnum,
     &      pointr,header)
            call pk_r8_opus('HOU',infovec(16),iend,mhl,errnum,
     &      pointr,header)
         endif
c DG090312
c TIDY UP!! - SIA/SIS for Wollongong, not within JFB's tlalevel system
         call pk_r8_opus('SIA',infovec(28),iend,mhl,errnum,
     &   pointr,header)
         call pk_r8_opus('SIS',infovec(29),
     &   iend,mhl,errnum,pointr,header)
         call pk_r8_opus('WSA',infovec(10),iend,mhl,errnum,
     &   pointr,header)
         call pk_r8_opus('WDA',infovec(13),iend,mhl,errnum,
     &   pointr,header)

c end DG090213

         if(tlalevel.ge.tla_ext_full) then
            call pk_r8_opus('WSA',infovec(10),iend,mhl,errnum,
     &      pointr,header)
            call pk_r8_opus('WSS',infovec(11),iend,mhl,errnum,
     &      pointr,header)
            call pk_r8_opus('WSM',infovec(12),iend,mhl,errnum,
     &      pointr,header)
            call pk_r8_opus('WDA',infovec(13),iend,mhl,errnum,
     &      pointr,header)
            call pk_r8_opus('WDS',infovec(14),iend,mhl,errnum,
     &      pointr,header)
            call pk_r8_opus('ZSA',infovec(17),iend,mhl,errnum,
     &      pointr,header)
            call pk_r8_opus('ZSS',infovec(18),iend,mhl,errnum,
     &      pointr,header)
            call pk_r8_opus('ZRM',infovec(20),iend,mhl,errnum,
     &      pointr,header)
            call pk_r8_opus('ZLM',infovec(21),iend,mhl,errnum,
     &      pointr,header)
            call pk_r8_opus('ZVA',infovec(22),iend,mhl,errnum,
     &      pointr,header)
            call pk_r8_opus('DAA',infovec(23),iend,mhl,errnum,
     &      pointr,header)
            call pk_r8_opus('DSM',infovec(24),iend,mhl,errnum,
     &      pointr,header)
            call pk_r8_opus('SAA',infovec(25),iend,mhl,errnum,
     &      pointr,header)
            call pk_r8_opus('SEA',infovec(26),iend,mhl,errnum,
     &      pointr,header)
            call pk_r8_opus('STM',infovec(27),iend,mhl,errnum,
     &      pointr,header)
            call pk_r8_opus('SIA',infovec(28),iend,mhl,errnum,
     &      pointr,header)
            call pk_r8_opus('SIS',infovec(29),iend,mhl,errnum,
     &      pointr,header)
            call pk_r8_opus('SOA',infovec(30),iend,mhl,errnum,
     &      pointr,header)
            call pk_r8_opus('SOE',infovec(31),iend,mhl,errnum,
     &      pointr,header)
            call pk_r8_opus('SDA',infovec(32),iend,mhl,errnum,
     &      pointr,header)
         endif
         call end_opus_prm(mhl,errnum,pointr,header)
         length=(pointr-bytstart)/4
         call put_i4_byte(dbborgpar,iend,mhl,errnum,dirpnt,header)
         call put_i4_byte(length,iend,mhl,errnum,dirpnt,header)
         call put_i4_byte(bytstart,iend,mhl,errnum,dirpnt,header)
         ndb=ndb+1
c
c  Create data status parameter block.  The block type variable 'idtags'
c  takes the following values:
c         1047 = Si spectrum
c         2071 = Si interferogram
c        33815 = InGaAs spectrum
c        34839 = InGaAs interferogram
c
         bytstart=pointr
         call pk_i4_opus('DPF',1,iend,mhl,errnum,pointr,header)
         call pk_i4_opus('NPT',maxind-minind+1,iend,mhl,errnum,
     &   pointr,header)
         call pk_r8_opus('FXV',fxv,iend,mhl,errnum,pointr,header)
         call pk_r8_opus('LXV',lxv,iend,mhl,errnum,pointr,header)
         call pk_r8_opus('CSF',1.d0,iend,mhl,errnum,pointr,header)
         call pk_r8_opus('MXY',mxy,iend,mhl,errnum,pointr,header)
         call pk_r8_opus('MNY',mny,iend,mhl,errnum,pointr,header)
         if(datype.eq.1) then
            call pk_st_opus('DXU','PNT',typ_enum,iend,
     &      mhl,errnum,pointr,header)
         elseif(datype.eq.2) then
            call pk_st_opus('DXU','WN',typ_enum,iend,
     &      mhl,errnum,pointr,header)
         endif
         write(stringa,'(2(i2.2,a1),i4)') day,'/',month,'/',year
         call pk_st_opus('DAT',stringa,typ_string,iend,
     &   mhl,errnum,pointr,header)
         write(stringa,'(3(i2.2,a1),i3.3,a8)')
     &   hour,':',minute,':',second,'.',millisec,' (UTC+0)'
         call pk_st_opus('TIM',stringa,typ_string,iend,
     &   mhl,errnum,pointr,header)
         call end_opus_prm(mhl,errnum,pointr,header)
         length=(pointr-bytstart)/4
         if(verbose.ge.4) then
            write(*,*)'idtags=',idtags
         endif
         call put_i4_byte(idtags,iend,mhl,errnum,dirpnt,header)
         call put_i4_byte(length,iend,mhl,errnum,dirpnt,header)
         call put_i4_byte(bytstart,iend,mhl,errnum,dirpnt,header)
         ndb=ndb+1
c
c  Prepare data block: byte reverse and write directory entry.  The block
c  type variable 'idtagd' takes the following values:
c         1031 = Si spectrum
c         2055 = Si interferogram
c        33799 = InGaAs spectrum
c        34823 = InGaAs interferogram
c
         if((iend.eq.bigendian).and.(errnum.eq.0)) then
            swapped=.true. 
            call rbyte(buffer(minind),4,maxind-minind+1)
         endif
         if(verbose.ge.4) then
            write(*,*)'idtagd=',idtagd
         endif
         call put_i4_byte(idtagd,iend,mhl,errnum,dirpnt,header)
         call put_i4_byte(maxind-minind+1,iend,mhl,errnum,dirpnt,header)
         call put_i4_byte(pointr,iend,mhl,errnum,dirpnt,header)
         ndb=ndb+1
c
c  Now that we've built all the blocks, we know the size of the
c  directory and we can put it back into the header block.
c
         call put_i4_byte(ndb,iend,mhl,errnum,ndbpnt,header)
      endif   ! (errnum.eq.0)
c
cc  If the file already exists, remove it.  Otherwise we would get an
cc  incorrect file length if the new one is shorter than the old one.
cc
c      lp=lnbc(path)
c      if(errnum.eq.0) then
c        inquire(file=path,exist=filexist,iostat=ierr)
c        if(ierr.ne.0) then
c          errnum=-10
c          write(*,'(2a)')'Error: inquire failed on file ',path(:lp)
c        elseif(filexist) then
c          if(verbose.ge.4) write(*,'(2a)')'Pre-deleting file ',path(:lp)
c          open(unit=luna,file=path,iostat=ierr)
c          if(ierr.ne.0) then
c            errnum=-10
c            write(*,'(2a)')'Error: Pre-delete failed on file ',path(:lp)
c          else
c            close(unit=luna,status='delete',iostat=ierr)
c            if(ierr.ne.0) then
c              errnum=-10
c            write(*,'(2a)')'Error: Pre-delete failed on file ',path(:lp)
c            endif
c          endif
c        endif
c      endif
c
c  Finally, write the header and data buffer in a single step.
c
c     write(*,*)'save_opus: pointr,minind,maxind=',pointr,minind,maxind

      lp=lnbc(path)
      if(errnum.eq.0) then
         open(unit=luna,file=path,form='unformatted',status='replace',
     &   access='direct',recl=pointr+4*(maxind-minind+1),iostat=ierr)
         if(ierr.ne.0) then
            errnum=-10
            write(*,'(2a)')'Error: open failed on file ',path(:lp)
         elseif(verbose.ge.3) then
            write(*,'(2a)')'Writing file: ',path(:lp)
            write(*,*)'....'
         endif
      endif
      if(errnum.eq.0) then
         write(luna,rec=1,iostat=ierr)
     &   (header(indexa),indexa=1,pointr),
     &   (buffer(indexa),indexa=minind,maxind)
         if(ierr.ne.0) then
            errnum=-10
            write(*,'(a)')'Error: file write failed'
         endif
      endif
      if(errnum.eq.0) then
         close(unit=luna,iostat=ierr)
         if(ierr.ne.0) then
            errnum=-10
            write(*,'(2a)')'Error: close failed on file ',path(:lp)
         endif
      endif
c  Undo the byte swap so that we don't modify the subroutine input argument
      if(swapped) call rbyte(buffer(minind),4,maxind-minind+1)

      return
      end
c===============================================================
      subroutine put_i4_byte(i4val,iendian,mhl,errnum,bytpoint,bythead)
c
c  Input:
c    i4val     I*4    Value of item to be stored in byte header
c    iendian   I*4    Endianness of computer
c    mhl       I*4    Maximum Header Length, in bytes
c
c  Input/Output:
c    errnum    I*4    Program error status
c    bytpoint  I*4    Pointer into byte header
c
c  Output:
c    bythead(mhl)BYTE Byte array to receive the item
c
      implicit none

      integer*4
     & i4val,      ! Subroutine input argument (see above)
     & iendian,    ! Subroutine input argument (see above)
     & mhl,        ! Subroutine input argument (see above)
     & errnum,     ! Subroutine input/output argument (see above)
     & bytpoint,   ! Subroutine input/output argument (see above)
     & bigendian,  ! Named constant for big endian detection 
     & i4loc,      ! Local copy of i4val
     & indexa      ! General loop index

      integer*1
     & bythead(mhl),!Subroutine output argument (see above)
     & bytloc(4)   ! Local copy of i4val as byte array

      parameter (bigendian=1)
      equivalence (i4loc,bytloc)

      if(errnum.ne.0) then
         continue
      elseif((bytpoint+4).ge.mhl) then
         errnum=-11
         write(*,'(a)')'Error in save_opus: increase parameter MHL'
      else
         i4loc=i4val
         if(iendian.eq.bigendian) call rbyte(i4loc,4,1)
         do indexa=1,4
            bytpoint=bytpoint+1
            bythead(bytpoint)=bytloc(indexa)
         enddo
      endif

      return
      end
c===============================================================
      subroutine put_r8_byte(r8val,iendian,mhl,errnum,bytpoint,bythead)
c
c  Input:
c    r8val     R*8    Value of item to be stored in byte header
c    iendian   I*4    Endianness of computer
c    mhl       I*4    Maximum Header Length, in bytes
c
c  Input/Output:
c    errnum    I*4    Program error status
c    bytpoint  I*4    Pointer into byte header
c
c  Output:
c    bythead(mhl)BYTE Byte array to receive the item
c
      implicit none

      integer*4
     & iendian,    ! Subroutine input argument (see above)
     & mhl,        ! Subroutine input argument (see above)
     & errnum,     ! Subroutine input/output argument (see above)
     & bytpoint,   ! Subroutine input/output argument (see above)
     & bigendian,  ! Named constant for big endian detection 
     & indexa      ! General loop index

      real*8
     & r8val,      ! Subroutine input argument (see above)
     & r8loc       ! Local copy of r8val

      integer*1
     & bythead(mhl),!Subroutine output argument (see above)
     & bytloc(8)   ! Local copy of r8val as byte array

      parameter (bigendian=1)
      equivalence (r8loc,bytloc)

      if(errnum.ne.0) then
         continue
      elseif((bytpoint+8).ge.mhl) then
         errnum=-11
         write(*,'(a)')'Error in save_opus: increase parameter MHL'
      else
         r8loc=r8val
         if(iendian.eq.bigendian) call rbyte(r8loc,8,1)
         do indexa=1,8
            bytpoint=bytpoint+1
            bythead(bytpoint)=bytloc(indexa)
         enddo
      endif

      return
      end
c===============================================================
      subroutine pk_i4_opus(param,i4val,iendian,mhl,errnum,
     & bytpoint,bythead)
c
c  Packs a 4-byte integer as an OPUS parameter into a byte header array.
c
c  Input:
c    param     C*3    Parameter name
c    i4val     I*4    Value of parameter
c    iendian   I*4    Endianness of computer
c    mhl       I*4    Maximum Header Length, in bytes
c
c  Input/Output:
c    errnum    I*4    Program error status
c    bytpoint  I*4    Pointer into byte header
c
c  Output:
c    bythead(mhl)BYTE Byte array to receive the item
c
      implicit none

      integer*4
     & i4val,      ! Subroutine input argument (see above)
     & iendian,    ! Subroutine input argument (see above)
     & mhl,        ! Subroutine input argument (see above)
     & errnum,     ! Subroutine input/output argument (see above)
     & bytpoint,   ! Subroutine input/output argument (see above)
     & bigendian,  ! Named constant for big endian detection 
     & i4loc,      ! Local copy of i4val
     & indexa      ! General loop index

      integer*2
     & i2loc(2),   ! Local variable to convert I*2 to bytes
     & typ,        ! Type of the OPUS parameter
     & rs          ! Reserved space for parameter, in 16-bit units

      integer*1
     & bythead(mhl),!Subroutine output argument (see above)
     & bytloc(4)   ! Local copy of i4val as byte array

      character
     & param*3     ! Subroutine input argument (see above)

      parameter (bigendian=1)
      parameter (typ=0, rs=2)
      equivalence (i4loc,i2loc,bytloc)

      if(errnum.ne.0) then
         continue
      elseif((bytpoint+12).ge.mhl) then
         errnum=-11
         write(*,'(a)')'Error in save_opus: increase parameter MHL'
      elseif(mod(bytpoint,2).ne.0) then
         errnum=-11
         write(*,'(a)')'Error in save_opus: pointer not on I*2 boundary'
      else
c
c  Save the null-terminated parameter name straight into the byte header.
c
         do indexa=1,3
            bytpoint=bytpoint+1
            call i4tobyte(ichar(param(indexa:indexa)),errnum,
     &      bythead(bytpoint))
         enddo
         bytpoint=bytpoint+1
         bythead(bytpoint)=0
c
c  Now handle the item type and reserved space.
c
         i2loc(1)=typ
         i2loc(2)=rs
         if(iendian.eq.bigendian) call rbyte(i2loc,2,2)
c         if(iendian.eq.bigendian) then
c            call i2rev(i2loc(1))
c            call i2rev(i2loc(2))
c         endif
         do indexa=1,4
            bytpoint=bytpoint+1
            bythead(bytpoint)=bytloc(indexa)
         enddo
c
c  Finally, save the parameter value.
c
         i4loc=i4val
         if(iendian.eq.bigendian) call rbyte(i4loc,4,1)
         do indexa=1,2*rs
            bytpoint=bytpoint+1
            bythead(bytpoint)=bytloc(indexa)
         enddo
      endif
 
      return
      end
c===============================================================
      subroutine pk_r8_opus(param,r8val,iendian,mhl,errnum,
     & bytpoint,bythead)
c
c  Packs an 8-byte real as an OPUS parameter into a byte header array.
c
c  Input:
c    param     C*3    Parameter name
c    r8val     I*4    Value of parameter
c    iendian   I*4    Endianness of computer
c    mhl       I*4    Maximum Header Length, in bytes
c
c  Input/Output:
c    errnum    I*4    Program error status
c    bytpoint  I*4    Pointer into byte header
c
c  Output:
c    bythead(mhl)BYTE Byte array to receive the item
c
      implicit none

      integer*4
     & iendian,    ! Subroutine input argument (see above)
     & mhl,        ! Subroutine input argument (see above)
     & errnum,     ! Subroutine input/output argument (see above)
     & bytpoint,   ! Subroutine input/output argument (see above)
     & bigendian,  ! Named constant for big endian detection 
     & indexa      ! General loop index

      integer*2
     & i2loc(2),   ! Local variable to convert I*2 to bytes
     & typ,        ! Type of the OPUS parameter
     & rs          ! Reserved space for parameter, in 16-bit units

      integer*1
     & bythead(mhl),!Subroutine output argument (see above)
     & bytloc(8)   ! Local copy of r8val as byte array

      real*8
     & r8val,      ! Subroutine input argument (see above)
     & r8loc       ! Local copy of r8val

      character
     & param*3     ! Subroutine input argument (see above)

      parameter (bigendian=1)
      parameter (typ=1, rs=4)
      equivalence (r8loc,i2loc,bytloc)

      if(errnum.ne.0) then
         continue
      elseif((bytpoint+16).ge.mhl) then
         errnum=-11
         write(*,'(a)')'Error in save_opus: increase parameter MHL'
      elseif(mod(bytpoint,2).ne.0) then
         errnum=-11
         write(*,'(a)')'Error in save_opus: pointer not on I*2 boundary'
      else
c
c  Save the null-terminated parameter name straight into the byte header.
c
         do indexa=1,3
            bytpoint=bytpoint+1
            call i4tobyte(ichar(param(indexa:indexa)),errnum,
     &      bythead(bytpoint))
         enddo
         bytpoint=bytpoint+1
         bythead(bytpoint)=0
c
c  Now handle the item type and reserved space.
c
         i2loc(1)=typ
         i2loc(2)=rs
         if(iendian.eq.bigendian) call rbyte(i2loc,2,2)
c           call i2rev(i2loc(1))
c           call i2rev(i2loc(2))
c         endif
         do indexa=1,4
            bytpoint=bytpoint+1
            bythead(bytpoint)=bytloc(indexa)
         enddo
c
c  Finally, save the parameter value.
c
         r8loc=r8val
         if(iendian.eq.bigendian) call rbyte(r8loc,8,1)
         do indexa=1,2*rs
            bytpoint=bytpoint+1
            bythead(bytpoint)=bytloc(indexa)
         enddo
      endif
 
      return
      end
c===============================================================
      subroutine pk_st_opus(param,cval,prmtyp,iendian,mhl,errnum,
     & bytpoint,bythead)
c
c  Packs a string as an OPUS parameter into a byte header array.
c
c  Input:
c    param     C*3    Parameter name
c    cval      C*(*)  Value of parameter
c    prmtyp    I*4    Parameter type (string/enum/senum)
c    iendian   I*4    Endianness of computer
c    mhl       I*4    Maximum Header Length, in bytes
c
c  Input/Output:
c    errnum    I*4    Program error status
c    bytpoint  I*4    Pointer into byte header
c
c  Output:
c    bythead(mhl)BYTE Byte array to receive the item
c
      implicit none

      integer*4
     & prmtyp,     ! Subroutine input argument (see above)
     & iendian,    ! Subroutine input argument (see above)
     & mhl,        ! Subroutine input argument (see above)
     & errnum,     ! Subroutine input/output argument (see above)
     & bytpoint,   ! Subroutine input/output argument (see above)
     & bigendian,  ! Named constant for big endian detection 
     & indexa,     ! General loop index
     & strlen,     ! String length
     & lnbc        ! Integer function Last Non-Blank Character in string

      integer*2
     & i2loc(2),   ! Local variable to convert I*2 to bytes
     & typ,        ! Type of the OPUS parameter
     & rs          ! Reserved space for parameter, in 16-bit units

      integer*1
     & bythead(mhl),!Subroutine output argument (see above)
     & bytloc(4)   ! Local copy of i2loc as byte array

      character
     & param*3,    ! Subroutine input argument (see above)
     & cval*(*)    ! Subroutine input argument (see above)

      parameter (bigendian=1)
      equivalence (i2loc,bytloc)

      if(errnum.ne.0) then
         continue
      elseif((bytpoint+12).ge.mhl) then
         errnum=-11
         write(*,'(a)')'Error in save_opus: increase parameter MHL'
      elseif(mod(bytpoint,2).ne.0) then
         errnum=-11
         write(*,'(a)')'Error in save_opus: pointer not on I*2 boundary'
      else
c
c  Compute length of input string.  The reserved space must take a
c  null terminator into account, hence the "+1".
c
         strlen=lnbc(cval)
         call i4tobyte(strlen/2+1,errnum,bytloc(1))  ! This silly thing to...
         rs=bytloc(1)                                ! avoid an ftnchek warning
c
c  Save the null-terminated parameter name straight into the byte header.
c
         do indexa=1,3
            bytpoint=bytpoint+1
            call i4tobyte(ichar(param(indexa:indexa)),errnum,
     &      bythead(bytpoint))
         enddo
         bytpoint=bytpoint+1
         bythead(bytpoint)=0
c
c  Now handle the item type and reserved space.
c
         call i4tobyte(prmtyp,errnum,bytloc(1))      ! This silly thing to...
         typ=bytloc(1)                               ! avoid an ftnchek warning
         i2loc(1)=typ
         i2loc(2)=rs
         if(iendian.eq.bigendian) call rbyte(i2loc,2,2)
c           call i2rev(i2loc(1))
c           call i2rev(i2loc(2))
c         endif
         do indexa=1,4
            bytpoint=bytpoint+1
            bythead(bytpoint)=bytloc(indexa)
         enddo
c
c  Save the parameter string value, along with a null terminator.
c
         do indexa=1,strlen
            bytpoint=bytpoint+1
            call i4tobyte(ichar(cval(indexa:indexa)),errnum,
     &      bythead(bytpoint))
         enddo
         bytpoint=bytpoint+1
         bythead(bytpoint)=0
c
c  Finally, make sure we end on an I*2 boundary
c
         if(mod(strlen,2).eq.0) then
            bytpoint=bytpoint+1
            bythead(bytpoint)=0
         endif
      endif
 
      return
      end
c===============================================================
      subroutine end_opus_prm(mhl,errnum,bytpoint,bythead)
c
c  Packs the 'END' mark of each parameter block into a byte header array.
c  Then ensure that the parameter block ends on a 4-byte boundary,
c  padding with zeros if necessary.
c
c  Input:
c    mhl       I*4    Maximum Header Length, in bytes
c
c  Input/Output:
c    errnum    I*4    Program error status
c    bytpoint  I*4    Pointer into byte header
c
c  Output:
c    bythead(mhl)BYTE Byte array to receive the item
c
      implicit none

      integer*4
     & mhl,        ! Subroutine input argument (see above)
     & errnum,     ! Subroutine input/output argument (see above)
     & bytpoint,   ! Subroutine input/output argument (see above)
     & indexa      ! General loop index

      integer*1
     & bythead(mhl)! Subroutine output argument (see above)

      if(errnum.ne.0) then
         continue
      elseif((bytpoint+10).ge.mhl) then
         errnum=-11
         write(*,'(a)')'Error in save_opus: increase parameter MHL'
      elseif(mod(bytpoint,2).ne.0) then
         errnum=-11
         write(*,'(a)')'Error in save_opus: pointer not on I*2 boundary'
      else
c
c  Save the 'END' straight into the byte header.
c
         bytpoint=bytpoint+1
         call i4tobyte(ichar('E'),errnum,bythead(bytpoint))
         bytpoint=bytpoint+1
         call i4tobyte(ichar('N'),errnum,bythead(bytpoint))
         bytpoint=bytpoint+1
         call i4tobyte(ichar('D'),errnum,bythead(bytpoint))
c
c  Now add the null termination plus the bytes normally occupied
c  by I*2 'type' and I*2 'reserved space': that's five bytes total.
c
         do indexa=1,5
            bytpoint=bytpoint+1
            bythead(bytpoint)=0
         enddo
c
c  If we are not on an I*4 boundary, add another two nulls.
c  Note that we were on an I*2 boundary on entry or we errored out.
c
         if(mod(bytpoint,4).ne.0) then
            bytpoint=bytpoint+1
            bythead(bytpoint)=0
            bytpoint=bytpoint+1
            bythead(bytpoint)=0
         endif
      endif
 
      return
      end
c===============================================================
      subroutine i4tobyte(i4var,errnum,bytevar)
c
c  After you stop rolling on the floor laughing, please email me
c  a better solution that still avoids the ftnchek warnings:
c  "intg expr XXXXX truncated to intg*1 YYYYY"    - JFB
c
      implicit none
      integer*4 i4var,errnum
      integer*1 bytevar,bytetabl(128)
      data bytetabl /0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,
     &               20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,
     &               36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,
     &               52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,
     &               68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,
     &               84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,
     &               100,101,102,103,104,105,106,107,108,109,110,111,
     &               112,113,114,115,116,117,118,119,120,121,122,123,
     &               124,125,126,127/

      if(errnum.ne.0) then
         continue
      elseif((i4var.le.127).and.(i4var.ge.0)) then
         bytevar=bytetabl(i4var+1)
      elseif(i4var.eq.176) then ! degree symbol - set to '
         bytevar=bytetabl(40)
      else
         errnum=-12
         write(*,'(a)')
     &   'Error: i4tobyte called with out-of-bound argument'
      endif

      return
      end
c===============================================================
      logical function is_leap_year(year)
      implicit none
      integer*4 year
      is_leap_year=((mod(year,4).eq.0).and.(mod(year,100).ne.0)).or.
     &              (mod(year,400).eq.0)
      return
      end
c===============================================================
      subroutine stardate_to_date (time,
     & year,month,day,hour,minute,second,millisec)
c
c  Convert from seconds since 1-Jan-2000 (stardate) to date and time.
c  NOTE: time before 1-Jan-2000 is not supported.
c
c  Input:
c    time      R*8    Time in seconds since 1-Jan-2000
c
c  Output:
c    year      I*4    Date component (full year, i.e. four digits, or more :) 
c    month     I*4    Date component
c    day       I*4    Date component
c    hour      I*4    Time component
c    minute    I*4    Time component
c    second    I*4    Time component
c    millisec  I*4    Time component
c
      implicit none

      real*8
     & time,       ! Subroutine input argument (see above)
     & timwrk      ! Temporary variable used in time computation

      integer*4
     & year,       ! Subroutine output argument (see above)
     & month,      ! Subroutine output argument (see above)
     & day,        ! Subroutine output argument (see above)
     & hour,       ! Subroutine output argument (see above)
     & minute,     ! Subroutine output argument (see above)
     & second,     ! Subroutine output argument (see above)
     & millisec,   ! Subroutine output argument (see above)
     & daysinyr,   ! Number of days in the current year
     & daysinmo    ! Number of days in the current month

      logical*4
     & is_leap_year,! Logical function
     & calculating  ! Loop control variable
c
c  Make a copy of the input value so that we don't modify it.
c
      timwrk=time
c
c  We first calculate the time, which is zero-based.
c
      millisec=nint((timwrk-dint(timwrk))*1000.d0)
      if(millisec.eq.1000) then
         millisec=0
         timwrk=timwrk+1.d0
      endif
      timwrk=dint(timwrk)/60.d0
      second=nint((timwrk-dint(timwrk))*60.d0)
      timwrk=dint(timwrk)/60.d0
      minute=nint((timwrk-dint(timwrk))*60.d0)
      timwrk=dint(timwrk)/24.d0
      hour=nint((timwrk-dint(timwrk))*24.d0)
c
c  Then we calculate the date, which is one-based, so we add 1.0 first.
c  At this point, 'timwrk' is in units of days.
c
      timwrk=dint(timwrk)+1.d0
c
c  Next we need to calculate the year because we must keep track of
c  the leap years and their effects on month and day of the month.
c  At the same time, we calculate the day number.  We're doing this
c  the "un-clever" way...  More elegant methods lead to embarrassing bugs.
c
      year=2000
      calculating=.true.
      do while(calculating)
         if(is_leap_year(year)) then
            daysinyr=366
         else
            daysinyr=365
         endif
         if(timwrk.gt.(dble(daysinyr)+0.5d0)) then
            year=year+1
            timwrk=timwrk-(dble(daysinyr))
         else
            calculating=.false.
         endif
      enddo
c
c  Finally we derive the month and day of the month, also being very
c  careful to avoid clever solutions.
c
      month=1
      calculating=.true.
      do while(calculating)
         if((month.eq.1).or.
     &      (month.eq.3).or.
     &      (month.eq.5).or.
     &      (month.eq.7).or.
     &      (month.eq.8).or.
     &      (month.eq.10).or.
     &      (month.eq.12)) then
            daysinmo=31
         elseif((month.eq.4).or.
     &      (month.eq.6).or.
     &      (month.eq.9).or.
     &      (month.eq.11)) then
            daysinmo=30
         elseif(is_leap_year(year)) then
            daysinmo=29
         else
            daysinmo=28
         endif
         if(timwrk.gt.(dble(daysinmo)+0.5d0)) then
            month=month+1
            timwrk=timwrk-(dble(daysinmo))
         else
            calculating=.false.
         endif
      enddo
      day=nint(timwrk)

      return
      end
