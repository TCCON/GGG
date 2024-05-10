      subroutine get_run_info(inpath,catyear,catmonth,catday,
     & catbatch,catslice,runno,verbose,mif,mns,nss,sfmcode,
     & stlimavg,stlimstd,runsta,timvec,infovec)

c  Input:
c    inpath      C*(*)  Directory path to info files
c    catyear     I*4    Year from catalog
c    catmonth    I*4    Month from catalog
c    catday      I*4    Day from catalog
c    catbatch    I*4    Batch (save set) number from catalog
c    catslice    I*4    Slice number from catalog, used in slice file name
c    runno       I*4    Run number increasing throughout the day
c    verbose     I*4    Level of verbosity for displayed messages
c    mif         I*4    Maximum number of info channels
c    mns         I*4    Maximum number of interferograms per scan set (max NSS)
c    nss         I*4    Number of Sample Scans in this set of slices
c    sfmcode     I*4    Numerical code to identify scan type (Sample ForM)
c    stlimavg    R*4    Limit for suntracker average intensity
c    stlimstd    R*4    Limit for suntracker intensity standard deviation

c  Input/output:
c    runsta(mns) I*4   Run starting slice number (if > 0) or run error (if < 0)
c    timvec(mns) R*8   Time vector contains one entry for each scan

c  Output:
c    infovec(mif,mns) R*8 Information produced by real-time algorithm
c
      implicit none

      integer*4 j,nlabel,idum,
     & catyear,    ! Subroutine input argument (see above)
     & catmonth,   ! Subroutine input argument (see above)
     & catday,     ! Subroutine input argument (see above)
     & catbatch,   ! Subroutine input argument (see above)
     & catslice,   ! Subroutine input argument (see above)
     & runno,      ! Subroutine input argument (see above)
     & verbose,    ! Subroutine input argument (see above)
     & mif,        ! Subroutine input argument (see above)
     & mns,        ! Subroutine input argument (see above)
     & nss,        ! Subroutine input argument (see above)
     & sfmcode,    ! Subroutine input argument (see above)
     & runsta(mns),! Subroutine input/output argument (see above)
     & luni,       ! Logical Unit Number for info file read
     & inpstat,    ! Value of IOSTAT from parameter input file read
     & errcode,    ! Error code (0=ok)
     & iscan,      ! Scan index
     & inlen,
     & indexa,     ! General loop index
     & fnbc,       ! Integer function First Non-Blank Character in string
     & lnbc        ! Integer function Last Non-Blank Character in string

      parameter (luni=22,nlabel=33)

      real*4
     & stlimavg,   ! Subroutine input argument (see above)
     & stlimstd    ! Subroutine input argument (see above)

      real*8
     & timvec(mns), ! Subroutine input/output argument (see above)
     & infovec(mif,mns),! Subroutine output argument (see above)
     & avgdrift     ! Average of time drifts

      character
     & label(nlabel)*22,
     & inpath*(*),  ! Subroutine input argument (see above)
     & dirsep*1,    ! Directory Seperator (/ or \)
     & filename*128,! Full file name for auxiliary info files
     & inpstr*999,  ! String used to read entries from input file
     & stringa*11,  ! String used to format integer display
     & stringb*11   ! String used to format integer/float display

      logical*4
     & filexist    ! Keeps track of file existence

      include 'header_indices.inc'

      data label/
     & "IFHum",
     & "IFSSrcT",
     & "IFS_P",
     & "ScBlkl_T",
     & "InGaAs_R",
     & "Si_R",
     & "Latitude",
     & "Longitude",
     & "Altitude",
     & "Zeno_WindSpeed_avg",
     & "Zeno_WindSpeed_std",
     & "Zeno_WindSpeed_max",
     & "Zeno_WindDir_avg",
     & "Zeno_WindDir_std",
     & "Zeno_Temp_avg",
     & "Zeno_RH_avg",
     & "Zeno_SolarRadiance_avg",
     & "Zeno_SolarRadiance_std",
     & "Zeno_Press_avg",
     & "Zeno_Rain_max",
     & "Zeno_Lightning_max",
     & "Zeno_VBatt_avg",
     & "Dome_azi_avg",
     & "Dome_Status_max",
     & "ST_tpg_azi_avg",
     & "ST_tpg_ele_avg",
     & "ST_TPS_max",
     & "ST_t_int_avg",
     & "ST_t_int_std",
     & "ST_off_azi_avg",
     & "ST_off_ele_avg",
     & "ST_Tdrift_avg",
     & "IFSDT_avg"/

      idum = src_off  ! Avoid compiler warning (unused)
      idum = src_sun  ! Avoid compiler warning (unused)
      idum = src_mir  ! Avoid compiler warning (unused)
      idum = src_nir  ! Avoid compiler warning (unused)
      idum = sfm_script  ! Avoid compiler warning (unused)

      idum=bms_caf2 ! Avoid compiler warnings (unused)
      idum=bms_kbr  ! Avoid compiler warnings (unused)
      idum=bms_quartz  ! Avoid compiler warnings (unused)
      idum=bms_sica  ! Avoid compiler warnings (unused)
      idum=i_aptval  ! Avoid compiler warnings (unused)
      idum=i_bbw   ! Avoid compiler warnings (unused)
      idum=i_bfw   ! Avoid compiler warnings (unused)
      idum=i_aqmcode  ! Avoid compiler warnings (unused)
      idum=i_bmscode  ! Avoid compiler warnings (unused)
      idum=i_dtccode  ! Avoid compiler warnings (unused)
      idum=i_inscode  ! Avoid compiler warnings (unused)
      idum=i_laserate ! Avoid compiler warnings (unused)

      idum=i_dur   ! Avoid compiler warnings (unused)
      idum=i_foc   ! Avoid compiler warnings (unused)
      idum=i_hpf   ! Avoid compiler warnings (unused)
      idum=i_gfw   ! Avoid compiler warnings (unused)
      idum=i_gbw   ! Avoid compiler warnings (unused)
      idum=i_lpf   ! Avoid compiler warnings (unused)
      idum=i_mvd   ! Avoid compiler warnings (unused)
      idum=i_p2a   ! Avoid compiler warnings (unused)
      idum=i_p2k   ! Avoid compiler warnings (unused)
      idum=i_p2l   ! Avoid compiler warnings (unused)
      idum=i_p2r   ! Avoid compiler warnings (unused)
      idum=i_pgn   ! Avoid compiler warnings (unused)
      idum=i_pka   ! Avoid compiler warnings (unused)
      idum=i_pra   ! Avoid compiler warnings (unused)
      idum=i_prl   ! Avoid compiler warnings (unused)
      idum=i_pkl   ! Avoid compiler warnings (unused)
      idum=i_rsn   ! Avoid compiler warnings (unused)
      idum=i_sgna  ! Avoid compiler warnings (unused)
      idum=i_sgnb  ! Avoid compiler warnings (unused)
      idum=i_srccode   ! Avoid compiler warnings (unused)
      idum=i_sfmcode   ! Avoid compiler warnings (unused)
      idum=i_ssm   ! Avoid compiler warnings (unused)
      idum=i_ssp   ! Avoid compiler warnings (unused)
      idum=i_zff   ! Avoid compiler warnings (unused)
      idum=i_lwn   ! Avoid compiler warnings (unused)
      idum=aqm_sd  ! Avoid compiler warnings (unused)
      idum=aqm_sf  ! Avoid compiler warnings (unused)
      idum=aqm_sn  ! Avoid compiler warnings (unused)

      inlen=lnbc(inpath)
      dirsep=inpath(inlen:inlen)

c
c  Fill critical elements of output vectors with "not found" values.
      do iscan=1,nss
         infovec(28,iscan)=-1.0d0       ! SIA  ST_t_int_avg  VALUE IMPORTANT
         infovec(29,iscan)=-1.0d0       ! SIS  ST_t_int_std  VALUE IMPORTANT
         infovec(33,iscan)=-9999999.0d0 ! IDA  IFSDT_avg IFS125_TIME - NTP_TIME
         infovec(35,iscan)=-1.0d0       ! ISS  ScanStatus    VALUE IMPORTANT
      enddo
c
c  Loop over the scans in this scan set to locate their info files.
c  We do this for scans that are not flagged in error by checking 'runsta'.
      do iscan=1,nss
         if(runsta(iscan).gt.0) then
            errcode=0
c            call build_info_name(inpath,catyear,catmonth,catday,
c     &      catbatch,catslice,iscan-1,filename)
            write(filename,'(a,3i2.2,a1,i0,a1,a4,2a1,i0,a1,i0,a5)')
     &      inpath(1:inlen),mod(catyear,100),catmonth,catday,'.',
     &      catbatch,dirsep,'scan',dirsep,'b',catslice,'.',iscan-1,
     &      '.info'
            inquire(file=filename,exist=filexist)
            if(filexist) then
               open(luni,file=filename,status='old')
               if(verbose.ge.3) then
                  write(*,'(a,i0,2a)') 'Run ',runno+iscan-1,
     &           ' auxiliary info: ',filename(1:lnbc(filename))
               endif
c
c  Reject runs that don't have info files
            else
               errcode=1
               runsta(iscan)=-105
c               call build_reject_label(catyear,catmonth,catday,
c     &         runno+iscan-1,rejlabl)
               write(*,'(a,i0,a,i4.4,i2.2,i2.2,a,a)') 'Reject: run ',
     &         runno+iscan-1,' of ',catyear,catmonth,catday,
     &         ' is missing info file ',filename(1:lnbc(filename))
            endif
c
c  If found, read and parse the contents of the info file.
            inpstat=0
c  While not at EOF and no errors...
            do while( inpstat.eq.0 .and. errcode.eq.0 )
               call read_input_line(luni,errcode,inpstat,inpstr)

               if(index(inpstr,"ScanType").gt.0) then
                  if(verbose.ge.4) write(*,*) inpstr(1:lnbc(inpstr))
                  if((index(inpstr,"Error").le.0).and.
     &            (((sfmcode.eq.sfm_solar).and.
     &            (index(inpstr,"Solar").le.0)).or.
     &            ((sfmcode.eq.sfm_cell).and.
     &            (index(inpstr,"Cell").le.0)).or.
     &            ((sfmcode.eq.sfm_idle).and.
     &            (index(inpstr,"IdleScan").le.0)).or.
     &            ((sfmcode.eq.sfm_aeros).and.
     &            (index(inpstr,"Aerosol").le.0))))then
                     if(verbose.ge.2) then
                        write(*,'(a,i0,2a)')
     &                  'Warning: scan type mismatch between codes ',
     &                  sfmcode,' and ',inpstr(1:lnbc(inpstr))
                     endif
                  endif
               endif

               if(index(inpstr,"ScanStatus").gt.0) then
                  if(verbose.ge.4) then
                     write(*,*) inpstr(1:lnbc(inpstr))
                  endif
                  if(index(inpstr,"OK").gt.0) then
                     infovec(35,iscan)=0.d0
                  elseif(index(inpstr,"Stress").gt.0) then
                     infovec(35,iscan)=1.d0
                  elseif(index(inpstr,"Error").gt.0) then
                     infovec(35,iscan)=2.d0
                  elseif(index(inpstr,"Reset").gt.0) then
                     infovec(35,iscan)=3.d0
                  endif
                  if((verbose.ge.2).and.
     &            (nint(infovec(35,iscan)).ne.0)) then
                     write(*,'(a)') 'Warning: scan status is not OK'
                  endif
               endif

c  Search through label(j) to find a match to inpstr. 
c  If found, read associated data value into infovec(j).
               do j=1,nlabel
                  if(index(inpstr,label(j)(:lnbc(label(j)))).gt.0) then
                     if(verbose.ge.4) write(*,*) inpstr(1:lnbc(inpstr))
                     read(inpstr(index(inpstr,":")+1:lnbc(inpstr)),*,
     &               iostat=inpstat) infovec(j,iscan)
                     if(inpstat.ne.0) then
                        if(verbose.ge.2) then
                           write(*,'(2a)')'Warning: fmt error:',label(j)
                        endif
                        errcode=2
                        inpstat=0
                     endif    !  if(inpstat.ne.0) then
                  endif    !  if(index(inpstr,label(j)).gt.0) then
               end do   !  do j=1,nlabel

c            if(index(inpstr,"IFHum").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(1,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &             'Warning: info file format error for IFHum'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"IFSSrcT").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(2,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &             'Warning: info file format error for IFSSrcT'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"IFS_P").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(3,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &             'Warning: info file format error for IFS_P'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"ScBlkl_T").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(4,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &             'Warning: info file format error for ScBlkl_T'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"InGaAs_R").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(5,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &             'Warning: info file format error for InGaAs_R'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"Si_R").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(6,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &             'Warning: info file format error for Si_R'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"Latitude").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(7,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &             'Warning: info file format error for Latitude'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"Longitude").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(8,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &             'Warning: info file format error for Longitude'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"Altitude").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(9,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &             'Warning: info file format error for Altitude'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"Zeno_WindSpeed_avg").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(10,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &          'Warning: info file format error for Zeno_WindSpeed_avg'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"Zeno_WindSpeed_std").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(11,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &          'Warning: info file format error for Zeno_WindSpeed_std'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"Zeno_WindSpeed_max").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(12,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &          'Warning: info file format error for Zeno_WindSpeed_max'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"Zeno_WindDir_avg").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(13,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &            'Warning: info file format error for Zeno_WindDir_avg'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"Zeno_WindDir_std").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(14,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &            'Warning: info file format error for Zeno_WindDir_std'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"Zeno_Temp_avg").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(15,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &             'Warning: info file format error for Zeno_Temp_avg'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"Zeno_RH_avg").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(16,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &             'Warning: info file format error for Zeno_RH_avg'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"Zeno_SolarRadiance_avg").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(17,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &      'Warning: info file format error for Zeno_SolarRadiance_avg'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"Zeno_SolarRadiance_std").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(18,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &      'Warning: info file format error for Zeno_SolarRadiance_std'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"Zeno_Press_avg").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(19,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &             'Warning: info file format error for Zeno_Press_avg'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"Zeno_Rain_max").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(20,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &             'Warning: info file format error for Zeno_Rain_max'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"Zeno_Lightning_max").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(21,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &          'Warning: info file format error for Zeno_Lightning_max'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"Zeno_VBatt_avg").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(22,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &             'Warning: info file format error for Zeno_VBatt_avg'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"Dome_azi_avg").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(23,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &             'Warning: info file format error for Dome_azi_avg'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"Dome_Status_max").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(24,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &             'Warning: info file format error for Dome_Status_max'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"ST_tpg_azi_avg").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(25,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &             'Warning: info file format error for ST_tpg_azi_avg'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"ST_tpg_ele_avg").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(26,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &             'Warning: info file format error for ST_tpg_ele_avg'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"ST_TPS_max").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(27,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &             'Warning: info file format error for ST_TPS_max'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"ST_t_int_avg").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(28,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &             'Warning: info file format error for ST_t_int_avg'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"ST_t_int_std").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(29,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &             'Warning: info file format error for ST_t_int_std'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"ST_off_azi_avg").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(30,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &             'Warning: info file format error for ST_off_azi_avg'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"ST_off_ele_avg").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(31,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &             'Warning: info file format error for ST_off_ele_avg'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"ST_Tdrift_avg").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(32,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &             'Warning: info file format error for ST_Tdrift_avg'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif
c
c            if(index(inpstr,"IFSDT_avg").gt.0) then
c              if(verbose.ge.4) then
c                write(*,*) inpstr(1:lnbc(inpstr))
c              endif
c              indexb=index(inpstr,":")
c              read(inpstr(indexb+1:lnbc(inpstr)),*,iostat=inpstat)
c     &         infovec(33,iscan)
c              if(inpstat.ne.0) then
c                if(verbose.ge.2) then
c                  write(*,'(a)')
c     &             'Warning: info file format error for IFSDT_avg'
c                endif
c                errcode=2
c                inpstat=0
c              endif
c            endif

            enddo   ! while((inpstat.eq.0).and. [...]
c
c  If the instrument is vented, replace the IFS125 pressure reading
c  with the more accurate Zeno value.  This is necessary for the
c  air-to-vacuum correction of the wavenumber scale.
            if((errcode.eq.0).and.(infovec(3,iscan).gt.400.d0)) then
               if(verbose.ge.2) then
                  write(*,'(2a)') 'Warning: instrument is vented, ',
     &            'replacing its pressure with external value'
               endif
               infovec(3,iscan)=infovec(19,iscan)
            endif
c
c  Reject scans based on solar tracker total intensity.
            if((sfmcode.eq.sfm_solar)
     &         .and.(nint(infovec(28,iscan)).ne.(-1))
     &         .and.(stlimavg.ne.0.0)
     &         .and.(infovec(28,iscan).le.dble(stlimavg))) then
               runsta(iscan)=-101
c              write(stringb,'(f11.1)') infovec(28,iscan)
c              write(stringc,'(f11.1)') stlimavg
c              call build_reject_label(catyear,catmonth,catday,
c     &        runno+iscan-1,rejlabl)
               write(*,'(a12,i0,a4,i4.4,2i2.2,a28,f11.1,a4,f11.1)')
     &        'Reject: run ',runno+iscan-1,' of ',catyear,catmonth,
     &         catday,' has solar intensity AVG = ',infovec(28,iscan),
     &        ' <= ',stlimavg
            endif

            if((sfmcode.eq.sfm_solar)
     &        .and.(nint(infovec(28,iscan)).ne.(-1))
     &        .and.(nint(infovec(29,iscan)).ne.(-1))
     &        .and.(stlimstd.ne.0.0)
     &        .and.(infovec(29,iscan).ge.
     &          (dble(stlimstd)*infovec(28,iscan)))) then
               runsta(iscan)=-102
c              write(stringb,'(f11.1)') infovec(29,iscan)
c              write(stringc,'(f11.1)') dble(stlimstd)*infovec(28,iscan)
c              call build_reject_label(catyear,catmonth,catday,
c     &        runno+iscan-1,rejlabl)
               write(*,'(a12,i0,a4,i4.4,2i2.2,a28,f11.1,a4,f11.1)')
     &        'Reject: run ',runno+iscan-1,' of ',catyear,catmonth,
     &         catday,' has solar intensity STD = ',infovec(29,iscan),
     &        ' >= ',dble(stlimstd)*infovec(28,iscan)
            endif

         endif
      enddo   ! iscan=1,nss
c
c  Check consistency of time drifts.
      do iscan=1,nss-1
         if((nint(infovec(33,iscan)).ne.(-9999999)).and.
     &      (nint(infovec(33,iscan+1)).ne.(-9999999)).and.
     &      (dabs(infovec(33,iscan)-infovec(33,iscan+1)).gt.2.d0))then
            write(stringa,'(f11.1)') infovec(33,iscan)
            write(stringb,'(f11.1)') infovec(33,iscan+1)
            if(dabs(infovec(33,iscan)).lt.
     &         dabs(infovec(33,iscan)+3600.d0)) then
               if(dabs(infovec(33,iscan)).lt.
     &         dabs(infovec(33,iscan+1))) then
                  infovec(33,iscan+1)=infovec(33,iscan)
               else
                  infovec(33,iscan)=infovec(33,iscan+1)
               endif
            else
               if(dabs(infovec(33,iscan)+3600.d0).lt.
     &         dabs(infovec(33,iscan+1)+3600.d0)) then
                  infovec(33,iscan+1)=infovec(33,iscan)
               else
                  infovec(33,iscan)=infovec(33,iscan+1)
               endif
            endif
            if(verbose.ge.2) then
               write(*,'(5a,i0)') 'Warning: inconsistent time drifts ',
     &         stringa(fnbc(stringa):),' and ',
     &         stringb(fnbc(stringb):),', using ',
     &         nint(infovec(33,iscan))
            endif
         endif
      enddo   ! iscan=1,nss-1
c
c  Compute average of individual scan time drifts.
      avgdrift=0.d0
      indexa=0
      do iscan=1,nss
         if(nint(infovec(33,iscan)).ne.(-9999999)) then
            avgdrift=avgdrift+infovec(33,iscan)
            indexa=indexa+1
         endif
      enddo
      if(indexa.ne.0) avgdrift=avgdrift/dble(indexa)
      if(verbose.ge.3) write(*,'(a,f7.1)')'Average time drift',avgdrift
c
c  Apply average time drift correction to the slices of the whole scan set.
      do indexa=1,nss
         timvec(indexa)=timvec(indexa)-avgdrift
      enddo

      return
      end
