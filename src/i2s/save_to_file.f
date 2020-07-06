      subroutine save_to_file(datype,flimit,filename,ichan,outfmt,
     & verbose,progver,mip,mif,mi4,mr8,buffer,counter,nside,
     & minfreq,maxfreq,delimit,minmax,nlong,nshort,time,
     & i4head,r8head,DTCstr,INSstr,sivcfreq,pco_leni,pco_threshi,
     & izpd,sivcflag,dclevel,fvsi_calc,zpa,frzpda,shbar,sherr,lsemode,
     & fpilha,infovec,tlalevel,fpsfname,errnum)
c
c  Input:
c    datype       I*4    Data type (1=interferogram, 2=spectrum)
c    flimit       C*(*)  Name of file containing the frequency limits
c    filename     C*(*)  Full file name for output files
c    ichan        I*4    Channel number (1=InGaAs=slave, 2=Si=master)
c    outfmt       I*4    Format selection for output file
c    verbose      I*4    Level of verbosity for displayed messages
c    progver      R*8    Program version (date) to be stored in OPUS header
c    mip          I*4    Maximum number of input points
c    mif          I*4    Maximum number of info channels
c    mi4          I*4    Maximum number of I*4 items in file header
c    mr8          I*4    Maximum number of R*8 items in file header
c    buffer(mip/2)  R*4    Buffer containing interferogram/spectrum data
c    counter      I*4    Count of points in buffer
c    minfreq      R*8    Frequency of the first data point in buffer
c    maxfreq      R*8    Frequency one spectral point past the last one in buffer
c    delimit      C*(*)  Field delimiter for ASCII output file
c    minmax       I*4    Count of min-max pairs in ASCII output file (0=no min-max)
c    nlong        I*4    Number of points in one FWD or REV scan
c    nshort       I*4    Number of points in one FWD or REV scan
c    time         R*8    Time in seconds since 1-Jan-2000
c    i4head(mi4)  I*4    Vector holding the I*4 header items
c    r8head(mr8)  R*8    Vector holding the R*8 header items
c    DTCstr       C*40   Variable to hold detector string
c    INSstr       C*40   Variable to hold instrument string
c    sivcfreq     R*8    SIV-Correction frequency (cm-1)
c    pco_threshi  R*8    SIV-Correction frequency (cm-1)
c    izpd         I*4    Point index of zero path difference
c    sivcflag     I*4    
c    dclevel      R*8    
c    zpa          R*8    ZPD interferogram amplitude (phase-corrected)
c    infovec(mif) R*8   Information produced by real-time algorithm
c    tlalevel     I*4    Level of non-Bruker header items (three-letter acronyms)
c
c  Input/Output:
c    errnum     I*4    Program error status
c
      implicit none

      integer*4 
     & nside,
     & idum,
     & datype,      ! Subroutine input argument (see above)
     & ichan,       ! Subroutine input argument (see above)
     & outfmt,      ! Subroutine input argument (see above)
     & verbose,     ! Subroutine input argument (see above)
     & mip,         ! Subroutine input argument (see above)
     & mif,         ! Subroutine input argument (see above)
     & mi4,         ! Subroutine input argument (see above)
     & mr8,         ! Subroutine input argument (see above)
     & counter,     ! Subroutine input argument (see above)
     & minmax,      ! Subroutine input argument (see above)
     & nlong,       ! Subroutine input argument (see above)
     & nshort,      ! Subroutine input argument (see above)
     & i4head(mi4), ! Subroutine input argument (see above)
     & izpd,        ! Subroutine input argument (see above)
     & sivcflag,    !
     & lsemode,     ! Laser sampling error type
     & pco_leni,    !
     & tlalevel,    ! Subroutine input argument (see above)
     & errnum,      ! Subroutine input/output argument (see above)
     & chanloc,     ! Local copy of 'channel' avoids ftnchek warnings
     & minind,      ! Point index of the lowest wavenumber to be saved
     & maxind       ! Point index of the highest wavenumber to be saved

      real*4
     & buffer(mip/2), ! Subroutine input argument (see above)
     & compute_snr,snr,
     & shbar,       ! Laser sampling error (LSE)
     & sherr,       ! Laser sampling error uncertainty (LSU)
     & fpilha       ! Fraction of power in lower half of alias

      real*8
     & dclevel,     ! DC interferogram signal level at ZPD
     & fvsi_calc,   ! FVSI calculated from the smoothed interferogram
     & zpa,         ! ZPD interferogram amplitude (phase-corrected)
     & frzpda,      ! ZPD interferogram dip (due to non-linearity)
     & progver,     ! Subroutine input argument (see above)
     & time,        ! Subroutine input argument (see above)
     & r8head(mr8), ! Subroutine input argument (see above)
     & minfreq,     ! Subroutine input argument (see above)
     & maxfreq,     ! Subroutine input argument (see above)
     & infovec(mif),!Subroutine input argument (see above)
     & minwav,      ! Lowest wavenumber to be saved in spectrum
     & maxwav,      ! Highest wavenumber to be saved in spectrum
     & sivcfreq,    ! SIV-Correction Frequency (cm-1)
     & pco_threshi, ! Phase Correction Operator: fractional intensity threshold
     & pointspac    ! Point spacing, in cm-1

      character
     & flimit*(*),  ! Subroutine input argument (see above)       
     & fpsfname*(*),
     & filename*(*),!Subroutine input argument (see above)       
     & delimit*(*), ! Subroutine input argument (see above)
     & DTCstr*(*),  ! String to hold detector description
     & INSstr*(*)   ! String to hold instrument description

      include 'header_indices.inc'

      idum=i_ssm        !  Avoid compiler warning (unused parameter)
      idum=i_ssp        !  Avoid compiler warning (unused parameter)
      idum=i_srccode    !  Avoid compiler warning (unused parameter)
      idum=i_sfmcode    !  Avoid compiler warning (unused parameter)
      idum=i_dtccode    !  Avoid compiler warning (unused parameter)
      idum=i_pgn        !  Avoid compiler warning (unused parameter)
      idum=i_sgna       !  Avoid compiler warning (unused parameter)
      idum=i_sgnb       !  Avoid compiler warning (unused parameter)
      idum=i_hpf        !  Avoid compiler warning (unused parameter)
      idum=i_gfw        !  Avoid compiler warning (unused parameter)
      idum=i_gbw        !  Avoid compiler warning (unused parameter)
      idum=i_bfw        !  Avoid compiler warning (unused parameter)
      idum=i_bbw        !  Avoid compiler warning (unused parameter)
      idum=i_rsn        !  Avoid compiler warning (unused parameter)
      idum=i_pka        !  Avoid compiler warning (unused parameter)
      idum=i_pkl        !  Avoid compiler warning (unused parameter)
      idum=i_pra        !  Avoid compiler warning (unused parameter)
      idum=i_prl        !  Avoid compiler warning (unused parameter)
      idum=i_p2a        !  Avoid compiler warning (unused parameter)
      idum=i_p2l        !  Avoid compiler warning (unused parameter)
      idum=i_p2r        !  Avoid compiler warning (unused parameter)
      idum=i_p2k        !  Avoid compiler warning (unused parameter)
      idum=i_zff        !  Avoid compiler warning (unused parameter)
      idum=i_inscode    !  Avoid compiler warning (unused parameter)
      idum=i_bmscode    !  Avoid compiler warning (unused parameter)
      idum=i_aqmcode    !  Avoid compiler warning (unused parameter)
      idum=i_laserate   !  Avoid compiler warning (unused parameter)
      idum=i_lwn        !  Avoid compiler warning (unused parameter)
      idum=i_foc        !  Avoid compiler warning (unused parameter)
      idum=i_aptval     !  Avoid compiler warning (unused parameter)
      idum=i_lpf        !  Avoid compiler warning (unused parameter)
      idum=i_dur        !  Avoid compiler warning (unused parameter)
      idum=i_mvd        !  Avoid compiler warning (unused parameter)
c GCT      idum=tla_ext_none !  Avoid compiler warning (unused parameter)
c GCT      idum=tla_ext_min  !  Avoid compiler warning (unused parameter)
c GCT      idum=tla_ext_pgr  !  Avoid compiler warning (unused parameter)
c GCT      idum=tla_ext_full !  Avoid compiler warning (unused parameter)
      idum=src_sun      !  Avoid compiler warning (unused parameter)
      idum=src_nir      !  Avoid compiler warning (unused parameter)
      idum=src_mir      !  Avoid compiler warning (unused parameter)
      idum=src_off      !  Avoid compiler warning (unused parameter)
      idum=sfm_solar    !  Avoid compiler warning (unused parameter)
      idum=sfm_cell     !  Avoid compiler warning (unused parameter)
      idum=sfm_idle     !  Avoid compiler warning (unused parameter)
      idum=sfm_script   !  Avoid compiler warning (unused parameter)
      idum=sfm_aeros    !  Avoid compiler warning (unused parameter)
      idum=bms_caf2     !  Avoid compiler warning (unused parameter)
      idum=bms_kbr      !  Avoid compiler warning (unused parameter)
      idum=bms_quartz   !  Avoid compiler warning (unused parameter)
      idum=bms_sica     !  Avoid compiler warning (unused parameter)
      idum=aqm_sd       !  Avoid compiler warning (unused parameter)
      idum=aqm_sf       !  Avoid compiler warning (unused parameter)
      idum=aqm_sn       !  Avoid compiler warning (unused parameter)
      idum=i_inscode    !  Avoid compiler warning (unused parameter)  

      chanloc=ichan      ! Avoids ftnchek warnings

      if(errnum.eq.0) then
         if(datype.eq.2) then  ! spectra

            snr=compute_snr(counter,buffer,16,0.25)
            write(*,*)'save_to_file: SNR=',snr

            call get_flimit(flimit,chanloc,errnum,minwav,maxwav)
            if(sivcfreq.gt.minwav) then
c              write(*,'(2a)')'Warning: SIV correction frequency',
c     &        ' exceeds starting frequency of spectrum'
c              write(*,*)sivcfreq,minwav
            endif

            if(errnum.eq.0) then
               pointspac=(maxfreq-minfreq)*nside/dble(counter)
               minind=nint((minwav-minfreq)/pointspac)+1
               if(minind.lt.1) minind=1
               maxind=nint((maxwav-minfreq)/pointspac)+1
               if(maxind.gt.counter/nside) maxind=counter/nside
               if(verbose.ge.4) then
                  write(*,*) 'minwav,maxwav,minind,maxind,counter=',
     &            minwav,maxwav,minind,maxind,counter
               endif
            endif
         else
            minind=1
            maxind=counter
         endif

         if(outfmt.eq.1) then
            call save_ascii(filename,verbose,mip/2,datype,buffer,
     &      counter,minfreq,maxfreq,minind,maxind,delimit,minmax,errnum)
         elseif(outfmt.eq.2) then
            call save_binary(filename,verbose,mip/2,buffer,
     &      minind,maxind,errnum)
         elseif(outfmt.eq.3) then
c           write(*,*)' save_opus: minwav, maxwav=',minwav,maxwav
c           write(*,*)' save_opus: minfreq, maxfreq=',minfreq,maxfreq
c           write(*,*)' save_opus: minind, maxind=',minind,maxind

            if (nside.eq.0) nside=1 
            call save_opus(filename,verbose,progver,mip/2,mif,mi4,mr8,
     &      datype,chanloc,buffer,counter/nside,
     &      minfreq,maxfreq,minind,maxind,
     &      nlong,nshort,time,i4head,r8head,DTCstr,INSstr,sivcfreq,
     &      pco_leni,pco_threshi,izpd,
     &      sivcflag,dclevel,fvsi_calc,zpa,frzpda,
     &      shbar,sherr,lsemode,fpilha,snr,
     &      infovec,tlalevel,fpsfname,errnum)
            write(*,*)' saved_opus: '
         endif    ! outfmt.eq.1
      endif
      return
      end
c ===================================================================
      subroutine get_flimit(flimit,ichan,errnum,minwav,maxwav)
c
c  Input:
c    flimit    C*(*)  Name of file containing the frequency limits
c    ichan     I*4    Channel number (1=InGaAs=slave, 2=Si=master)
c
c  Input/Output:
c    errnum    I*4    Program error status
c
c  Output:
c    minwav    R*8    Lowest wavenumber to be saved in spectrum
c    maxwav    R*8    Highest wavenumber to be saved in spectrum
c
      implicit none

      integer*4
     & ichan,      ! Subroutine input argument (see above)
     & errnum,     ! Subroutine input/output argument (see above)
     & inpstat,    ! Value of IOSTAT from parameter input file read
     & ijk,     ! General loop index
     & lnbc,       ! Integer function Last Non-Blank Character in string
     & luna        ! Logical Unit Number for flimit file

      real*8
     & minwav,     ! Subroutine output argument (see above)
     & maxwav      ! Subroutine output argument (see above)

      character
     & flimit*(*), ! Subroutine input argument (see above)       
     & inpstr*999  ! String used to read entries from input file

      parameter (luna=19)

      if(errnum.eq.0) then
         open(unit=luna,file=flimit,status='old',iostat=inpstat)
         if(inpstat.ne.0) then
            errnum=-9
            write(*,'(2a)')'Error: open failed on file ',
     &      flimit(1:lnbc(flimit))
         endif
         ijk=0
         do while((inpstat.eq.0).and.(errnum.eq.0).and.
     &   (ijk.lt.ichan))
            call read_input_line(luna,errnum,inpstat,inpstr)
            if((errnum.eq.0).and.(inpstat.ne.0)) then
               write(*,'(a)') 'Error in flimit file: too few entries'
               errnum=-9
            elseif(errnum.eq.0) then
               read(inpstr,*,iostat=inpstat) minwav,maxwav
               if(inpstat.ne.0) then
                  write(*,'(a)')'Error in flimit file: bad entry format'
                  errnum=-9
               elseif(maxwav.lt.minwav) then
                  write(*,'(a)') 'Error in flimit file: bad ordering'
                  errnum=-9
               else
                  ijk=ijk+1
               endif
            endif
         enddo
         close(unit=luna,iostat=inpstat)
         if(inpstat.ne.0) then
            errnum=-9
            write(*,'(2a)')'Error: close failed on file ',
     &      flimit(1:lnbc(flimit))
         endif
      endif

      return
      end
c ===================================================================
      subroutine save_ascii(path,verbose,mip,datype,buffer,counter,
     & minfreq,maxfreq,minind,maxind,delimit,minmax,errnum)
c
c  Input:
c    path      C*(*)  Path to ASCII output file
c    verbose   I*4    Level of verbosity for displayed messages
c    mip       I*4    Maximum number of input points
c    datype    I*4    Data type (1=interferogram, 2=spectrum)
c    buffer(mip)R*4   Buffer containing interferogram/spectrum data
c    counter   I*4    Count of points in buffer
c    minfreq   R*8    Frequency of the first data point in buffer
c    maxfreq   R*8    Frequency one spectral point past the last one in buffer
c    minind    I*4    Point index of the lowest wavenumber to be saved
c    maxind    I*4    Point index of the highest wavenumber to be saved
c    delimit   C*(*)  Field delimiter for ASCII output file
c    minmax    I*4    Count of min-max pairs in ASCII output file (0=no min-max)
c
c  Input/Output:
c    errnum    I*4    Program error status
c
      implicit none

      integer*4
     & verbose,    ! Subroutine input argument (see above)
     & mip,        ! Subroutine input argument (see above)
     & datype,     ! Subroutine input argument (see above)
     & counter,    ! Subroutine input argument (see above)
     & minind,     ! Subroutine input argument (see above)
     & maxind,     ! Subroutine input argument (see above)
     & minmax,     ! Subroutine input argument (see above)
     & errnum,     ! Subroutine input/output argument (see above)
     & luna,       ! Logical Unit Number for ASCII output file
     & ierr,       ! Error status returned by file I/O
     & fnbc,       ! Integer function First Non-Blank Character in string
     & lnbc,       ! Integer function Last Non-Blank Character in string
     & ijk,     ! General loop index
     & mmcnt,      ! Number of input points per min/max pair
     & first       ! Starting point index for the current min/max pair

      real*4
     & buffer(mip),! Subroutine input argument (see above)
     & rmin,       ! Minimum Y-value for the current min/max pair
     & rmax        ! Maximum Y-value for the current min/max pair

      real*8
     & minfreq,    ! Subroutine input argument (see above)
     & maxfreq,    ! Subroutine input argument (see above)
     & pointspac   ! Point spacing, in cm-1

      character
     & path*(*),   ! Subroutine input argument (see above)
     & delimit*(*),! Subroutine input argument (see above)
     & strx*13, ! String used to format X-value display
     & stry*13  ! String used to format Y-value display

      parameter (luna=19)

      if(errnum.eq.0) then
         open(unit=luna,file=path,iostat=ierr)
         if(ierr.ne.0) then
            errnum=-10
            write(*,'(2a)')'Error: open failed on file ',
     &      path(1:lnbc(path))
         elseif(verbose.ge.3) then
            write(*,'(2a)')'Writing file: ',path(1:lnbc(path))
         endif
      endif
      if(errnum.eq.0) then
         pointspac=(maxfreq-minfreq)/dble(counter)
         if((minmax.eq.0).and.(datype.eq.1)) then
            do ijk=minind,maxind
               if(errnum.eq.0) then
                  write(strx,'(i13)') ijk
                  write(stry,'(f10.7)') buffer(ijk)
                  write(luna,'(3a)',iostat=ierr) strx(fnbc(strx):),
     &            delimit(1:1),stry(fnbc(stry):)
                  if(ierr.ne.0) then
                     errnum=-10
                     write(*,'(a)')'Error: file write failed'
                  endif
               endif
            enddo  ! ijk=minind,maxind
         elseif(minmax.eq.0) then
            do ijk=minind,maxind
               if(errnum.eq.0) then
                  write(strx,'(f13.8)')minfreq+(dble(ijk-1)*pointspac)
                  write(stry,'(f13.9)') buffer(ijk)
                  write(luna,'(3a)',iostat=ierr) strx(fnbc(strx):),
     &            delimit(1:1),stry(fnbc(stry):)
                  if(ierr.ne.0) then
                     errnum=-10
                     write(*,'(a)')'Error: file write failed'
                  endif
               endif
            enddo  ! ijk=minind,maxind
         else
            mmcnt=(maxind-minind)/minmax
            ijk=minind
            do while(ijk.le.maxind)
               first=ijk
               rmin=buffer(ijk)
               rmax=rmin
               ijk=ijk+1
               do while((ijk.le.maxind).and.
     &         ((ijk-first).lt.mmcnt))
                  if(buffer(ijk).lt.rmin) rmin=buffer(ijk)
                  if(buffer(ijk).gt.rmax) rmax=buffer(ijk)
                  ijk=ijk+1
               enddo ! while((ijk.le.maxind).and.((ijk-....
               if(datype.eq.1) then
                  if(errnum.eq.0) then
                     write(strx,'(i13)') (first+ijk)/2
                     write(stry,'(f10.7)') rmin
                     write(luna,'(3a)',iostat=ierr) strx(fnbc(strx):),
     &               delimit(1:1),stry(fnbc(stry):)
                     if(ierr.ne.0) then
                        errnum=-10
                        write(*,'(a)')'Error: file write failed'
                     endif
                  endif
                  if(errnum.eq.0) then
                     write(stry,'(f10.7)') rmax
                     write(luna,'(3a)',iostat=ierr) strx(fnbc(strx):),
     &               delimit(1:1),stry(fnbc(stry):)
                     if(ierr.ne.0) then
                        errnum=-10
                        write(*,'(a)')'Error: file write failed'
                     endif
                  endif
               else
                  if(errnum.eq.0) then
                     write(strx,'(f13.6)')
     &               minfreq+(dble((first+ijk-2)/2)*pointspac)
                     write(stry,'(f13.4)') rmin
                     write(luna,'(3a)',iostat=ierr) strx(fnbc(strx):),
     &               delimit(1:1),stry(fnbc(stry):)
                     if(ierr.ne.0) then
                        errnum=-10
                        write(*,'(a)')'Error: file write failed'
                     endif
                  endif
                  if(errnum.eq.0) then
                     write(stry,'(f13.4)') rmax
                     write(luna,'(3a)',iostat=ierr) strx(fnbc(strx):),
     &               delimit(1:1),stry(fnbc(stry):)
                     if(ierr.ne.0) then
                        errnum=-10
                        write(*,'(a)')'Error: file write failed'
                     endif
                  endif
               endif
            enddo  ! while(ijk.le.maxind)
         endif
         close(unit=luna,iostat=ierr)
         if(ierr.ne.0) then
            errnum=-10
            write(*,'(2a)')'Error: close failed on file ',
     &      path(1:lnbc(path))
         endif
      endif
      return
      end
c ===================================================================
      subroutine save_binary(path,verbose,mip,buffer,
     & minind,maxind,errnum)
c
c  Input:
c    path      C*(*)  Path to binary output file
c    verbose   I*4    Level of verbosity for displayed messages
c    mip       I*4    Maximum number of input points
c    buffer(mip)R*4   Buffer containing interferogram/spectrum data
c    minind    I*4    Point index of the lowest wavenumber to be saved
c    maxind    I*4    Point index of the highest wavenumber to be saved
c
c  Input/Output:
c    errnum    I*4    Program error status
c
      implicit none

      integer*4
     & lp,
     & verbose,    ! Subroutine input argument (see above)
     & mip,        ! Subroutine input argument (see above)
     & minind,     ! Subroutine input argument (see above)
     & maxind,     ! Subroutine input argument (see above)
     & errnum,     ! Subroutine input/output argument (see above)
     & iendian,    ! Endianness of computer
     & bigendian,  ! Named constant for big endian detection 
     & luna,       ! Logical Unit Number for binary output file
     & ierr,       ! Error status returned by file I/O
     & lnbc,       ! Integer function Last Non-Blank Character in string
     & ijk      ! General loop index

      parameter (bigendian = 1)

      real*4
     & buffer(mip) ! Subroutine input argument (see above)

      character
     & path*(*)    ! Subroutine input argument (see above)

      logical*4
c     & filexist,   ! Keeps track of file existence
     & swapped     ! Indicates if we byte-swapped the data vector

      parameter (luna=19)
c
c  Initialize variables.
      swapped=.false. 
c
c  Byte-swap the input vector on big-endian machines.
      if(errnum.eq.0) then
         call getendian(iendian)
         if(iendian.eq.bigendian) then
            swapped=.true. 
            call rbyte(buffer(minind),4,maxind-minind+1)
         endif
      endif
c
cc  If the file already exists, remove it.  Otherwise we would get an
cc  incorrect file length if the new one is shorter than the old one.
c        lp=lnbc(path)
c        inquire(file=path,exist=filexist,iostat=ierr)
c        if(ierr.ne.0) then
c          errnum=-10
c          write(*,'(2a)')'Error: inquire failed on file ',path(:lp)
c        elseif(filexist) then
c          if(verbose.ge.4) write(*,'(2a)')'Pre-deleting file ',path(:lp)
c          open(unit=luna,file=path,iostat=ierr)
c          if(ierr.ne.0) then
c            errnum=-10
c            write(*,'(2a)')'Error: Pre-delete failed  file ',path(:lp)
c          else
c            close(unit=luna,status='delete',iostat=ierr)
c            if(ierr.ne.0) then
c              errnum=-10
c              write(*,'(2a)')'Error: Pre-delete failed, file ',path(:lp)
c            endif
c          endif
c        endif
c
c  Now write the file.
      lp=lnbc(path)
      if(errnum.eq.0) then
         open(unit=luna,file=path,access='direct',status='replace',
     &   form='unformatted',recl=4*(maxind-minind+1),iostat=ierr)
         if(ierr.ne.0) then
            errnum=-10
            write(*,'(2a)')'Error: open failed on file ',
     &      path(1:lnbc(path))
         elseif(verbose.ge.3) then
            write(*,'(2a)')'Writing file: ',path(1:lnbc(path))
         endif
      endif
      if(errnum.eq.0) then
         write(luna,rec=1,iostat=ierr) (buffer(ijk),ijk=minind,maxind)
         if(ierr.ne.0) then
            errnum=-10
            write(*,'(a)')'Error: file write failed'
         endif
      endif
      if(errnum.eq.0) then
         close(unit=luna,iostat=ierr)
         if(ierr.ne.0) then
            errnum=-10
            write(*,'(2a)')'Error: close failed on file ',
     &      path(1:lnbc(path))
         endif
      endif
c
c  Undo the byte swap so that we don't modify the subroutine input argument.
      if(swapped) call rbyte(buffer(minind),4,maxind-minind+1)

      return
      end
