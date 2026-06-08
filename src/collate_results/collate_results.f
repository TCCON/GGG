c  Program collate_results
c
c  Reads the .col output files (gas_1234.runlog.col) produced by GFIT,
c  collates them with auxiliary data from the runlog, and then writes
c  the information to a spreadsheet-format, XYPLOT-readable, output file
c  named  runlog.tsw, runlog.vsw, etc. 
c
c  INPUT FILES:
c     multiggg.sh          batch file containing names of .col files
c     runlog.xrl           runlog file
c     gas_1234.runlog.col  files containing column amounts
c
c  OUTPUT FILES:
c     runlog.xsw           Spreadsheet of individual window values
c                          for each spectrum

c  Also does a weighted average of the rms residual
c     RMS_bar =  Sum(wt.rmsocl)/Sum(wt)
c  wt = 1/(rmsocl^2+0.1)
c
c  If rmsocl is very large, then wt is small and that measurement
c  is essentially ignored. If rmsocl is very small, which is common
c  in completely saturated windows in which GFIT has increased CL,
c  that measurement has no more weight than a measurement with a
c  moderate value of rmsocl.
c
c  Runlog is opened, read, and closed for every window, avoiding
c  having to store in memory, apart from the spectrum names.
c
c  One record of the output file is written for each runlog record
c  having a different time or name.  Runlog records with the same
c  time/names are combined. Missing values are written if a .col
c  file doesn't cover particular spectra, maintaining synchrony
c  with windows that do.
c
c  Info in the .col file that you want in the auxiliary part of the
c  output file and that does not appear in runlog (e.g. zmin) must
c  be read into a vector.

c  Pseudo-code:
c      do icol=1,ncol     !  main loop over windows/col files
c         open(lunr_col,file=colfile,status='old') ! .col file
c         open(lunr_rl,file=rlgfile, status='old')   !DG000906
c         call read_runlog_header(lunr_rl,data_fmt_read_rl,col_labels_rl)
c         iobs=1
c         do while (iobs.lt.mobs)  ! Loop over RL spectra with different times or names
c            read(lunr_col,'(a)',end=84) col_string
c            read(col_string,csformat) specname_col,zpdtim,..vsf,vsf_err
c84          continue
c            do while(specname_col.ne.specname_rlg)
c               call read_runlog_data_record(lunr_rl,....
c               if(istat.ne.0) go to 14     ! Exit Loop  iobs=1,mobs
c               if(col1.eq.':' .or. col1.eq.';') cycle
c               delta_t=zpdtim-zpdwas
c
cc  Create separate entries in the .vsw file for spectra whose
cc  ZPD times differ by more than MAX_DELTA_T
c               if(iyr.ne.iyrwas .or. doy.ne.doywas .or.
c     &            dabs(delta_t).ge.max_delta_t ) then  ! New observation
c
cc  Force separate entries in the .xsw file if names are different.
c                  if(specname_rlg(4:lg).ne.specname_was(4:lg)) then
c                     c1sign(iobs)=col1
c                     spectrum(iobs)=specname_rlg
c                     year(iobs)=r8year 
c                     yaux(icol,iobs)=
c                     jval=jval+ncol    !  jval=icol+ncol*(iobs-1)
c                     yobs(jval)=ymiss  !  initialize to missing value
c                     yerr(jval)=ymiss  !  initialize to missing value
c                     iobs=iobs+1
c                  endif  !   if(specname_rlg.ne.specname_was
c               else
cc                  write(*,*)'not new obs '//specname_rlg(:20),delta_t
c               endif   ! if(iyr.ne.iyrwas .or. doy.ne.doywas 
c               iyrwas=iyr
c               doywas=doy
c               zpdwas=zpdtim
c               specname_was=specname_rlg
c            end do     ! while(specname_col.ne.specname_rlg)
c            c1sign(iobs)=col1
c            spectrum(iobs)=specname_rlg
c            year(iobs)=r8year 
c            yaux(icol,iobs)=
c            yobs(jval)=
c            yerr(jval)=
c            ncount(iobs)=ncount(iobs)+1
c            npp=npp+1
c         end do   !  while(iobs.lt.mobs) !  Loop over spectra
c         stop 'Increase parameter MROW'
c14       nobs=iobs  ! the number of records in the .col and .xsw files
c         close(lunr_rl)
c         close(lunr_col)
c      end do     ! icol=1,ncol    main loop over windows/.col files
c
c  Write out all analyzed abundances to the .?sw disk file.
c      open(lunw_xsw,file=outfile,status='unknown')
c      do iobs=1,nobs
c         write(lunw_xsw,....(yobs(),yerr()
c      end do  !  iobs=1,row
c      close(lunw_xsw)

      implicit none
      include "../gfit/ggg_int_params.f"
      include "params.f"
      include "../comn/postproc_params.f"

      integer iobs,jtg,j,k,ktg,kt0,lnbc,lloc,idot,lr,
     & iggg,igog,lsf,sf_flag,
     & npp,ldot,totnit,ntc,mlabel,mval,jval,lunw_nts,
     & lunr_mul,lunr_rl,lunr_col,lunw_rpt,lunw_xsw,lunw_sites,
     & lunw_trp,mcol,ncol,icol,mobs,lnit,lg,ld,tnt,tnd,
     & msite,nsite,jsite,lunw_miss,lunw_vsf,lunw_fssg,idum,iunder,
     & nauxcol,nlhead,nmiss,naddn,iyrwas,doywas,mask1,mask2,
     & nobs,nit,mit,i,iyr,doy,istat,nfound
      parameter (lunw_sites=21)    ! collate_sites.rpt
      parameter (lunr_mul=51)      ! multiggg.sh 
      parameter (lunr_rl=53)       ! runlog
      parameter (lunr_col=55)      ! .col file
      parameter (lunw_xsw=58)      ! output file (.xsw)
      parameter (lunw_rpt=60)      ! .rpt file
      parameter (lunw_trp=61)      ! .trp file (transpose of .xsw)
      parameter (lunw_nts=62)      ! .nts file (negative time step)
      parameter (lunw_miss=65)     ! .missing file
      parameter (lunw_vsf=71)      ! .vsf file
      parameter (lunw_fssg=72)     ! .fssg file
      parameter (mcol=628)         ! Total number of columns/windows
      parameter (mobs=360000)      ! Max number of output observations
      parameter (mval=35100000)    ! Max number of values (NROW*NCOL)
      parameter (mlabel=17200)     ! Max # of characters in col labels
      parameter (msite=20000)      ! Max number of observation sites

c  nobs is the number of observations = records in the .xsw output file
c  For single-detector acquisition, this is the number of records in RLG.
c  For dual-detector acquisition, this is half the number of records in RLG.

      character 
     & dl*1,
     & col1*1,
     & apf*2,
     & rlgfile*(mfilepath),
     & version*64,
     & data_fmt_read_rl*256,col_labels_rl*320,
     & ans*1, mode*4,
     & collabel*(mlabel),
     & sflabel*(mlabel),
     & auxcol*(13*mauxcol),
     & outfile*80,
     & specname_rlg*(nchar),
     & c1sign(mobs)*1,
     & spectrum(mobs)*(nchar),
     & output_fmt*40,
     & windows(mcol)*12,
     & gggdir(mpath),
     & specname_was*(nchar),
     & chardum*16,
     & addn_lines(maddln)*(mcharhead)

      real*8 asza,
     & zobs,
     & rmin,rmax,vmin,vmax,vemax,vemin,tvsf,tvse,ww,
     & graw,obslat,obslon,opd,
     & r8was,r8year,r8ydiff,year(mobs),
     & sitelat(msite),
     & sitelon(msite),
     & sitealt(msite),
     & wt,twt,trms,lasf,wavtkr,aipl,sia,fvsi,azim,wspd,wdir,osds

      real*4
     & fcen(mcol),width(mcol),zmin(mobs),app,
     & airmass,rmsoclpc,fs,r2,sg,sfval,
     & vsf,vsf_err,ymiss,small,cl,tilt,cc,zlo,ovcol,
c     & rversion,
     & yaux(mauxcol,mobs),
     & yobs(mval), yerr(mval), qc(mobs),ytot(mobs)
      parameter (ymiss=9.8765e+35,small=1.0E-18)

      real*8 zpdtim,tout,pout,hout,tins,pins,hins,fovi,fovo,amal,
     & snr,zenoff,zoff,zpdwas,max_delta_t,delta_t

      integer*4 bytepw,ifirst,ilast,possp,ntot(msite),nday(msite),
     & ncount(mobs),prepend_specname

      character yaux_string*640,command_string*256

      logical is_em27
      logical append_qcflag
      logical isclose_s ! function to test value equivalence within float erro

      data ntot/msite*0/
      data nday/msite*0/
      data ncount/mobs*0/

      idum=mcolvav  ! Prevent compiler warning (unused parameter)
      idum=mgas     ! Prevent compiler warning (unused parameter)
      idum=mlev     ! Prevent compiler warning (unused parameter)
c      idum=mobs_qc  ! Prevent compiler warning (unused parameter)
      idum=mspeci   ! Prevent compiler warning (unused parameter)
      idum=mvmode   ! Prevent compiler warning (unused parameter)
      idum=ncell    ! Prevent compiler warning (unused parameter)
      idum=mcolvsw  ! Prevent compiler warning (unused parameter)
      idum=mrow_qc  ! Prevent compiler warning (unused parameter)
      chardum=countfmt ! Prevent compiler warning (unused parameter)

      data zmin/mobs*0.0/
      nsite=0 ! initialize nsite
      sf_flag=0
      append_qcflag=.false.
      prepend_specname=0  ! prevent compiler warning

      call get_ggg_environment(gggdir,dl)
      version=
     &' collate_results          Version 2.10    2020-07-31   GCT,JLL'
      write(6,*) version
      lr=0

c NB: Geoff's version adds "poff" after "solzen". Will probably need
c to remerge that at some point. -JLL
      auxcol=
     &'      year         day          hour         run    '//
     &'      lat          long         zobs         zmin   '//
     &'      solzen       azim         osds         opd    '//
     &'      fovi         amal         graw         tins   '//
     &'      pins         tout         pout         hout   '//
     &'      sia          fvsi         wspd         wdir   '
c     &'     year          day           hour          run           lat 
c     &         long          zobs          zmin          asza          
c     &poff          azim          osds          opd           fovi     
c     &   amal          graw          tins          pins          tout 
c     &      pout          hout          sia           fvsi          wspd
c     &        wdir'
      call substr(auxcol,cdum,1,nauxcol)

c  Initialize character arrays (Necessary for the G77 compiler).
      do i=1,mlabel
         collabel(i:i)=' '
         sflabel(i:i)=' '
      end do

      if (iargc() == 0) then
         write(6,'(a)')
     $' vsf [t], original col (o), vertical col [v], los col [l],',
     $' continuum [c], tilt [m], cc [n],  freq shift [f],', 
     $'  solar-gas freq shift [s], rms [r]: '
         read(5,*) ans
      elseif(iargc() .ge. 1) then
         call getarg(1, ans)
      else
         stop 'Use: $gggpath/bin/collate_results t/v/l/c/m/n/f/s/r'
      endif

      is_em27 = .false.
      if (iargc() .ge. 2) then
         call getarg(2, mode)
         if (mode .eq. 'em27') then
            is_em27 = .true.
            write(*,*) 'Using EM27 settings'
         endif
      endif

c  Find the number of windows/columns (NCOL)
      open(lunr_mul,file='multiggg.sh',status='old')
      do icol=1,mcol  
121      read(lunr_mul,'(a)',end=99) tabel
         if(tabel(1:1).eq.':') go to 121
         if(lnbc(tabel).le.0) go to 121
      end do  ! icol=1,mcol     !  loop (over windows)
      read(lunr_mul,*,end=99) tabel
      write(*,*) 'mcol=',mcol
      stop 'Increase parameter mcol'
 99   close(lunr_mul)
      ncol=icol-1

      open(lunw_vsf,file='collate_results.vsf',status='unknown')
      write(lunw_vsf,*) 2, 6
      write(lunw_vsf,*)'       icol        ispec     freq            wid
     &             vsf             vsf_error'
      open(lunw_fssg,file='collate_results.fssg',status='unknown')
      write(lunw_fssg,*) 2, 9
      write(lunw_fssg,*)' icol  ispec       year        freq            
     & wid         graw         osds           fs           sg'
c  Read in the retrieved absorber amounts (YOBS+-YERR)
      open(lunr_mul,file='multiggg.sh',status='old')
      open(lunw_rpt,file='collate_results.rpt',status='unknown')
      open(lunw_nts,file='collate_results.nts',status='unknown')
      write(lunw_rpt,*)2,16
      write(lunw_rpt,*) 'iwin  fcen  fcen_error  Nrow    Npp   NIT   '//
     &'%Conv  RMS_Min  RMS_Mean  RMS_Max   VSF_min    VSF_bar   '//
     &'VSF_bar_error  VSF_max     VERR_min    VERR_max'
      do icol=1,ncol     !  main loop (over windows / .col files)
         npp=0
         twt=small
         trms=0.0d0
         totnit=0
         ntc=0
         jval=icol-ncol
135      read(lunr_mul,'(a)') tabel
         if(tabel(1:1).eq.':') go to 135
         if(lnbc(tabel).le.0)  go to 135
         iggg=index(tabel,'.ggg')
         igog=index(tabel,'.gog')
         colfile=tabel(index(tabel,' ')+1:max(iggg,igog))//'col'
         write(*,*) 'Opening '//colfile(:64)
         open(lunr_col,file=colfile,status='old') ! .col file
         idot=index(colfile,'.')
         collabel=collabel(:lnbc(collabel)+2)//colfile(:idot-1)//' '//
     $   colfile(:idot-1)//'_error'
         windows(icol)=colfile(:idot-1)
c         write(*,*) 'icol,windows=',icol,windows(icol)
         iunder=index(colfile(:idot-1),'_')
         if(iunder.gt.0) idot=iunder
c
c  Read header lines of .col file and locate column containing
c  "OVC_gas" in order to read data from appropriate target gas.
         read(lunr_col,*) nlhead
         read(lunr_col,'(a)') gfit_version
         read(lunr_col,'(a)') gsetup_version
         do k=4,nlhead-2
            read(lunr_col,'(34x,a)')header_string
            if(k.eq.6) rlgfile=header_string(:150)    ! GCT 2009-03-04
            if(index(header_string,'runlog').gt.0)
     &      rlgfile=header_string(:mpath)
         end do
         csformat=header_string(:lnbc(header_string))
c         write(*,'(a)') 'csformat='//csformat

         read(lunr_col,'(a)') command_string
         lsf=index(command_string,'sf=')
         if(lsf.gt.0) then
            sf_flag=1
            read(command_string(lsf+3:),*) sfval
         else
            sfval=1.0
         endif
         write(sflabel(1+8*(icol-1):),'(f8.3)')sfval

c         write(*,'(a)') 'command_string='//command_string
         read(command_string,*) fcen(icol), width(icol), mit
c         write(*,*) 'icol,fcen(icol) = ',icol,fcen(icol),width(icol),
c     &  colfile(:idot-1) 
         read(lunr_col,'(a)')header_string
c         write(*,'(a)') 'header_string='//header_string
         lnit= index(header_string,'Nit')
c         write(*,*) 'lnit= ',lnit

         kt0=index(header_string,' OVC_')
         if(colfile(1:1).eq.'m'.or.colfile(1:1).eq.'v') then ! Kludge for InSb+InGaAs/InGaAs+Si
           ktg=1+index(header_string,'OVC_'//colfile(2:idot-1))
         else
           ktg=1+index(header_string,' OVC_'//colfile(:idot-1))
         endif
c         write(*,*) 'ktg= ',ktg,' OVC_'//colfile(:idot-1)
         ktg=1+(ktg-kt0)/42
c        ktg=1+(ktg-80)/48
c         write(*,*) 'ktg= ',ktg,' OVC_'//colfile(:idot-1)
c         write(*,*)'colfile,ktg=',colfile(:idot-1),ktg
         iyrwas=-99999
         doywas=-99999
         zpdwas=-99999.9d0
         rmin=999.9999
         rmax=-99.9999
         vmin=+1.E+38
         vmax=-1.E+38
         vemin=+1.E+38
         vemax=-1.E+38
         tvsf=small
         tvse=small
         specname_rlg=' '
         specname_was='x'
         nsite=0

         lr=lnbc(rlgfile)
         if (is_em27) then
            max_delta_t=0.00014 ! 0.5 s, recommended by J. Hedelius for EM27s
         elseif(rlgfile(lr-2:lr-2).eq.'o') then
            max_delta_t=0.0004  ! 1.44s (ACE)
         else                 ! Ground-based
c            max_delta_t=0.0025  ! 9.0s 
            max_delta_t=0.0014  ! 5.0s  GCT 2018-08-07
         endif

c  Add spectrum name to output files only on non-MkIV gnd data
         ld=lloc(rlgfile,dl)
         if(rlgfile(lr-2:lr-2).eq.'g' .and.
     &      rlgfile(ld+1:ld+3).ne.'m4_' .and.
     &      rlgfile(ld+1:ld+9).ne.'solar_all' .and.
     &      rlgfile(ld+1:ld+5).ne.'synth' .and.
     &      rlgfile(ld+1:ld+10).ne.'ll20101005' .and.
     &      rlgfile(ld+1:ld+4).ne.'mkiv' ) then
            prepend_specname = 1
         endif
c
c  Read auxilliary measurements from runlog
c         write(*,*) 'runlog=',rlgfile
         open(lunr_rl,file=rlgfile, status='old')   !DG000906
         call read_runlog_header(lunr_rl,data_fmt_read_rl,col_labels_rl)
         r8was=-9999999.9d0
         iobs=1
         do while (iobs.lt.mobs)  !  Loop over observation times in RLG
            specname_col='='
            read(lunr_col,'(a)',end=84) col_string
            if (lnbc(col_string).le.2) goto 84 ! skip blank line at EOF

            read(col_string,csformat)
     &      specname_col(:lnit-3),nit,cl,tilt,cc,fs,
     &      sg,zlo,rmsoclpc,zmin(iobs),
     &      (airmass,ovcol,vsf,vsf_err,jtg=1,ktg)

c            write(*,*)'icol,iobs=',icol,iobs,specname_col(:24)
            totnit=totnit+nit
            if(nit.lt.mit) ntc=ntc+1  ! Number of Times Converged
            if(rmsoclpc.le.0.0) then
               write(*,*) 'rmsoclpc <= 0  ',colfile,iobs,lnit
c               rmsoclpc=0.0001
               stop 'rmsoclpc <= 0'   ! Commented 2009-03-18
            endif
            if(rmsoclpc.gt.rmax) rmax=rmsoclpc
            if(rmsoclpc.lt.rmin) rmin=rmsoclpc
            if(vsf.lt.vmin) vmin=vsf
            if(vsf.gt.vmax) vmax=vsf
            if(vsf_err.lt.vemin) vemin=vsf_err
            if(vsf_err.gt.vemax) vemax=vsf_err
            ww=1.0d0/(0.000001+vsf_err**2)
            tvse=tvse+ww
            tvsf=tvsf+vsf*ww
84          continue

c            write(*,*) icol,iobs,'specname_col.ne.specname_rlg: |'
c     &     //specname_col(:24)//'|'//specname_rlg(:24)//'|'
            do while(specname_col.ne.specname_rlg)
               call read_runlog_data_record(lunr_rl,data_fmt_read_rl,
     &         col1,specname_rlg,iyr,doy,zpdtim,obslat,obslon,zobs,
     &         asza,zenoff,azim,osds,opd,fovi,fovo,amal,ifirst,ilast,
     &         graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     &         tout,pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)
c               write(*,*)'read_runlog: iobs,istat=',iobs,istat,
c     &         specname_rlg
c               if(istat.ne.0)write(*,*)'Called readrunlog: istat=',istat
               if(istat.ne.0) go to 14     ! Exit Loop  iobs=1,mobs
               if(col1.eq.':' .or. col1.eq.';') cycle

               do jsite=1,nsite
                  if(abs(obslat-sitelat(jsite)) .lt. small .and.
     &               abs(obslon-sitelon(jsite)) .lt. small .and.
     &               abs(zobs-sitealt(jsite))   .lt. small) go to 66
               end do

c         New site
               nsite=nsite+1
               sitelat(jsite)=obslat
               sitelon(jsite)=obslon
               sitealt(jsite)=zobs

66             ntot(jsite)=ntot(jsite)+1
               if(iyr.ne.iyrwas .or. doy.ne.doywas) then
                  nday(jsite)=nday(jsite)+1
c                  if(jsite.eq.9)write(56,*)iobs,zobs,iyr,doy,nday(jsite)
               endif

               lg=lnbc(specname_rlg)
               ldot=index(specname_rlg,'.')
               delta_t=zpdtim-zpdwas
c               write(*,*)lg,ldot,specname_rlg(:lg),delta_t, max_delta_t
c
c  Create separate entries in the .vsw file for spectra whose
c  ZPD times differ by more than MAX_DELTA_T
c               write(*,*)'iyrwas,iyr=',iyrwas,iyr
c               write(*,*)'doywas,doy=',doywas,doy
               if(iyr.ne.iyrwas .or. doy.ne.doywas .or.
     &            dabs(delta_t).ge.max_delta_t ) then  ! New observation
c
c  The following if-statement shouldn't be necessary. But occasionally
c  you get simultaneous InGaAs/Si scans with very different ZPD times.
c  You only want them to have separate entries in the .vsw file if their
c  names are different. So mask out the parts of the spectrum name that
c  are allowed to differ. If everything else is the same, merge into
c  single entry.
                  if(ldot.eq.17) then ! TCCON spectrum naming convention
                     mask1=ldot-1
                     mask2=ldot-1
                  else if(ldot.eq.9) then ! MkIV spectrum naming convention
                     mask1=2
                     mask2=3
                  else
c Other spectrum naming convention; this will just check to see whether
c the spectrum name is identical (i.e., a do-nothing check)
                     mask1=2
                     mask2=1
                  end if
c                  if(specname_rlg(4:lg).ne.specname_was(4:lg)
c                  if(specname_rlg(4:ldot-2).ne.specname_was(4:ldot-2)
c     &           .or. specname_rlg(ldot:).ne.specname_was(ldot:)
                  if(specname_rlg(:mask1-1)//specname_rlg(mask2+1:lg)
     &            .ne.
     &            specname_was(:mask1-1)//specname_was(mask2+1:lg)) then
                     c1sign(iobs)=col1
                     spectrum(iobs)=specname_rlg
c  Code works because 2000 was leap year (it won't work for 2100)
                     if(mod(iyr,4).eq.0) then  ! Leap Year
                        r8year=iyr+(doy+zpdtim/24.0d0)/366.0d0
                     else                     ! non-Leap_year
                        r8year=iyr+(doy+zpdtim/24.0d0)/365.0d0
                     endif
                     year(iobs)=r8year ! added by DW 20120105
c JLL - removed this instance of writing the aux variables because it
c was causing zmin to be written incorrectly for the first spectrum
c of an InGaAs-only block in a mixed InGaAs-InSb runlog.
c                     write(yaux_string,*)doy+zpdtim/24.d0,zpdtim,
c     &               iobs,obslat,obslon,zobs,zmin(iobs),asza+zenoff,
c     &               azim,osds,opd,fovi,amal,graw,tins,pins,
c     &               tout,pout,hout,sia,fvsi,wspd,wdir
c                     read(yaux_string,*) (yaux(j,iobs),j=2,nauxcol)
                     jval=jval+ncol   !  jval=icol+ncol*(iobs-1)
                     yobs(jval)=ymiss
                     yerr(jval)=ymiss
                     iobs=iobs+1
                  endif  !   if(specname_rlg.ne.specname_was
               else
c                  write(*,*)'not new obs '//specname_rlg(:20),delta_t
               endif   ! if(iyr.ne.iyrwas .or. doy.ne.doywas 
               iyrwas=iyr
               doywas=doy
               zpdwas=zpdtim
               specname_was=specname_rlg
            end do     ! while(specname_col.ne.specname_rlg)
c            write(lunw_vsf,*)icol,iobs-1,fcen(icol),width(icol)/2,
c     &      vsf,vsf_err
        
c            iobs=iobs+1
            c1sign(iobs-1)=col1
            spectrum(iobs-1)=specname_rlg
            ncount(iobs-1)=ncount(iobs)+1
c            write(*,*) icol,iobs,ncount(iobs),spectrum(iobs)
c            r8year=iyr+(doy+zpdtim/24.0d0)/366.0d0
c  Code works because 2000 was leap year (it won't work for 2100)
            if(mod(iyr,4).eq.0) then  ! Leap Year
               r8year=iyr+(doy+zpdtim/24.0d0)/366.0d0
            else                     ! non-Leap_year
               r8year=iyr+(doy+zpdtim/24.0d0)/365.0d0
            endif
            r8ydiff=r8year-r8was
c  Report negative time-steps in the runlogs times.
c  Need only do this for the first window (avoid repetition).
            if( r8ydiff .lt. -0.00000001d0 .and. icol.eq.1) then
               write(lunw_nts,'(a,a,2f12.6)')
     &      '  Negative time step (runlog unsorted?) ',
     &         specname_rlg,r8was,r8year
            endif
            r8was=r8year
            year(iobs-1)=r8year

            write(lunw_fssg,'(2i6,7f14.7)')icol,iobs-1,year(iobs-1),
     &      fcen(icol),width(icol)/2,graw,osds,fs,sg

c NB: to return to being consistent with Geoff code, we will need to
c write asza and zenoff separately in the future. For now, I am keeping
c them added together so that the GGG2020 .xsw format does not change.
c -JLL
            write(yaux_string,*)doy+zpdtim/24.d0,zpdtim,iobs-1,
     &      obslat,obslon,zobs,zmin(iobs-1),asza+zenoff,azim,osds,opd,
     &      fovi,amal,graw,tins,pins,tout,pout,
     &      hout,sia,fvsi,wspd,wdir
            read(yaux_string,*) (yaux(j,iobs-1),j=2,nauxcol)
c               write(*,*)'2:asza,zenoff=',asza,zenoff,asza+zenoff
c               write(*,*)'yaux(9)=',yaux(9,iobs)

            if(ans.eq.'t') then
               yobs(jval)=vsf
               yerr(jval)=vsf_err
            elseif(ans.eq.'o') then
               yobs(jval)=ovcol
               yerr(jval)=ovcol*1.e-6
            elseif(ans.eq.'v') then
               yobs(jval)=vsf*ovcol
               yerr(jval)=vsf_err*ovcol
            elseif(ans.eq.'l') then
               yobs(jval)=vsf*ovcol*airmass
               yerr(jval)=vsf_err*ovcol*airmass
            elseif(ans.eq.'f') then
               yobs(jval)=fs
               yerr(jval)=vsf_err
c              yerr(jval)=vsf_err*abs(yobs(icol,iobs))
            elseif(ans.eq.'s') then
               yobs(jval)=sg
               yerr(jval)=vsf_err
            elseif(ans.eq.'m') then
               yobs(jval)=tilt
               yerr(jval)=vsf_err
            elseif(ans.eq.'n') then
               yobs(jval)=cc
               yerr(jval)=vsf_err
            elseif(ans.eq.'r') then
               yobs(jval)=rmsoclpc
               yerr(jval)=vsf_err*rmsoclpc
            elseif(ans.eq.'c') then
               yobs(jval)=cl
               yerr(jval)=cl*rmsoclpc
            else
               stop 'unknown option'
            endif
            if(yerr(jval).gt.3.4028E+38) yerr(jval)=3.4028E+38
            if(yerr(jval).lt.1.175E-38) yerr(jval)=1.175E-38
            npp=npp+1
            r2=3.0*abs(rmsoclpc)  ! rmsoclpc typically ~ 0.33%
c            wt= 1./(r2+1.0/r2)
            wt= 1./(r2+0.1)
            twt=twt+wt
            trms=trms+wt*rmsoclpc
         end do   !  iobs.lt.mobs  !  Loop over observation times in RLG
         stop 'Increase parameter MOBS'
14       nobs=iobs-1  ! the number of records in the .xsw file.
         if(ncol*nobs.gt.mval) then
            write(*,*)ncol,nobs
            write(*,*)'Increase parameter MVAL to ',nobs*ncol
            stop 'collate'
         endif
         close(lunr_rl)
         close(lunr_col)
         app=npp+small
         write(lunw_rpt,'(i3,f9.2,f8.2,2i8,2f7.2,3f9.4,6(1pe12.4))')
     &   icol,fcen(icol),width(icol)/2,npp,nobs,float(totnit)/app,
     &   100*float(ntc)/app,rmin,trms/twt,rmax,vmin,tvsf/tvse,
     &   sqrt(tvse*nobs)/tvsf,vmax,vemin,vemax
      end do  ! icol=1,mcol     !  main loop (over windows)
      close(lunr_mul)
      close(lunw_rpt)
      close(lunw_nts)
      close(lunw_vsf)
      close(lunw_fssg)
c====================================================================
c      do k=lr,1,-1
c         if(ichar(rlgfile(k:k)) .eq. 92) go to 101  ! backslash
c         if(ichar(rlgfile(k:k)) .eq. 47) go to 101  ! forward slash
c      end do
c101   outfile=rlgfile(k+1:lr-3)//ans//'sw'
      outfile=rlgfile(lloc(rlgfile,dl)+1:lr-3)//ans//'sw'
c==================================================================
      if (append_qcflag) then
         auxcol=auxcol(:lnbc(auxcol))//'      qcflag '
         nauxcol=nauxcol+1
         call generate_qc_flag(nobs,spectrum,qc)
         do iobs=1, nobs
            yaux(nauxcol,iobs)=qc(iobs)
         enddo
      endif
c====================================================================
c      call substr(auxcol,cdum,1,nauxcol)
c      if(nauxcol.gt.mauxcol) then
c          write(*,*)' mauxcol, nauxcol = ',mauxcol,nauxcol
c          stop 'increase parameter mauxcol'
c      endif

      output_fmt='(a1,f13.8,NNf13.5,MMMM(1pe13.5))'
      write(output_fmt(11:12),'(i2.2)') nauxcol-1
      write(output_fmt(19:22),'(i4.4)') 2*ncol
      if (prepend_specname.ge.1) then
         output_fmt='(a57,'//output_fmt(2:)
         auxcol='  spectrum   '//auxcol(:13*nauxcol)
         nauxcol=nauxcol+1
      endif
      write(*,*) ' Output format = '//output_fmt

c  Write out all analyzed abundances to the .?sw disk file.
      open(lunw_xsw,file=outfile,status='unknown')
      addn_lines(1) = version(:lnbc(version))
      addn_lines(2) = gfit_version(:lnbc(gfit_version))
      addn_lines(3) = gsetup_version(:lnbc(gsetup_version))
      naddn = 3

c  If there are defined scale factors for our windows, then
c  we need to add those into the header. Add them before the
c  missing and format lines so that those two are always the
c  last two lines for read_postproc_header.
      if(sf_flag .eq. 1) then
        addn_lines(4) = 'sf='//sflabel(:lnbc(sflabel))
        naddn = 4
      endif
      call write_postproc_header(lunw_xsw, nauxcol+2*ncol, nobs, 
     & nauxcol, dble(ymiss), output_fmt, addn_lines, naddn, 0)
      write(lunw_xsw,'(a)') auxcol(:lnbc(auxcol))//'  '//
     &   collabel(:lnbc(collabel)+1)

c  Note that the SIGN array and the following IF statement are merely
c  to support both the new and the old runlog formats.
      open(lunw_miss,file='collate_results.missing',status='unknown')
      nfound=0
      nmiss=0
      jval=0
      do iobs=1,nobs
         do k=1,ncol
            jval=jval+1     !    jval=k+ncol*(iobs-1)
c            if( yobs(jval).eq.ymiss .and.
c     &          yerr(jval).eq.ymiss) then
            if( isclose_s(yobs(jval), ymiss) .and.
     &          isclose_s(yerr(jval), ymiss) ) then
               nmiss=nmiss+1
               write(lunw_miss,*)'Missing: ',windows(k),
     &       '  '//spectrum(iobs),yobs(jval),yerr(jval)
            else
               nfound=nfound+1
            endif
         end do
         if (prepend_specname.ge.1) then
            write(lunw_xsw,output_fmt) spectrum(iobs),c1sign(iobs),
     &      year(iobs),(yaux(k,iobs),k=2,nauxcol-1),
     &      (yobs(k+ncol*(iobs-1)),yerr(k+ncol*(iobs-1)),k=1,ncol)
         else
            write(lunw_xsw,output_fmt) c1sign(iobs),year(iobs),
     &      (yaux(k,iobs),k=2,nauxcol),
     &      (yobs(k+ncol*(iobs-1)),yerr(k+ncol*(iobs-1)),k=1,ncol)
         endif

      end do  !  iobs=1,nobs
      close(lunw_miss)
      close(lunw_xsw)
c====================================================================
      do j=1,nobs
         ytot(j)=0.0
      end do
c  Write another output file containing the transpose of data matrix
      open(lunw_trp,file=outfile(:lnbc(outfile))//'.transpose',
     & status='unknown')
      write(lunw_trp,*)4,nobs+3
      write(lunw_trp,*)'Missing:',ymiss
      write(lunw_trp,'(a24,600000a36)')' Window fcen fcen_error ',
     & (spectrum(j),j=1,nobs)
      write(lunw_trp,'(a24,600000i4)')' Window fcen fcen_error ',
     & (j,j=1,nobs)
      do k=1,ncol
         write(lunw_trp,*)k,fcen(k),1+width(k)/2,
     &  (yobs(k+ncol*(j-1)),j=1,nobs)
         do j=1,nobs
           if(ytot(j).lt.ymiss) ytot(j)=ytot(j)+yobs(k+ncol*(j-1))
         end do
      end do
      close(lunw_trp)

c      write(67,*) 2,3
c      write(67,*) ' j zmin rms'
c      do j=1,nobs
c         write(67,*) j, yaux(8,j), ytot(j)/ncol
c      end do
c====================================================================
      write(6,*)
      write(6,*) outfile(:lnbc(outfile))//' contains:'
      write(6,'(i8,a)') prepend_specname,' spectrum name column'
      write(6,'(i8,a)') nauxcol,' auxiliary columns (nauxcol)'
      write(6,'(i8,a)') ncol,' data column pairs (value + error)'
      write(6,'(i8,a)') prepend_specname+nauxcol+2*ncol,' total columns'
      write(6,'(i8,a)') nobs,' data rows'
      write(6,'(i8,a,f5.1,a1)') nfound,' found values   = ',
     & 100*float(nfound)/(nfound+nmiss),'%'
      write(6,'(i8,a,f5.1,a30)') nmiss,' missing values = ',
     & 100*float(nmiss)/(nfound+nmiss),'%  See collate_results.missing'
c====================================================================
c  Report site information
      tnt=0
      tnd=0
      open(lunw_sites,file='collate_sites.rpt',status='unknown')
      write(lunw_sites,*)2,6
      write(lunw_sites,*)' #     Nobs    Nday     Lat     Long      Alt'
      do jsite=1,nsite
         write(lunw_sites,'(i3,2i8,3f12.3)') jsite,ntot(jsite)/ncol,
     &   nday(jsite)/ncol,sitelat(jsite),sitelon(jsite),sitealt(jsite)
         tnt=tnt+ntot(jsite)
         tnd=tnd+nday(jsite)
      end do
      write(lunw_sites,*)'Ntot=',tnt/ncol
      write(lunw_sites,*)'Nday=',tnd/ncol
      close(lunw_sites)
      stop
      end

c      subroutine set_yaux(y(2),y(3),y(4),y(5),y(6),y(7),y(8),y(9),y(10),
c     &y(11),y(12),y(13),y(14),y(15),y(16),y(17),y(18),y(19),y(20),y(21),
c     &y(22),y(23),y)
c      real*4 y(23)
c      return
c      end
