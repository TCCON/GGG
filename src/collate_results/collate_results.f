c  Program collate_results
c
c  Reads the .col output files (gas_1234.runlog.col) produced by GFIT,
c  collates them with auxiliary data from the runlog, and then writes
c  the information to a spreadsheet-format XYPLOT-readable output file
c  named  runlog.tsw, runlog.csw, etc.
c
c  INPUT FILES:
c     multiggg.sh          batch file containing names of .col files
c     runlog.xrl           runlog file
c     gas_1234.runlog.col  files containing column amounts
c
c  OUTPUT FILES:
c     runlog.xsw           Spreadsheet of individual window values
c                          for each spectrum

c  Does a weighted average of the rms residual
c     RMS_bar =  Sum(wt.rmsocl)/Sum(wt)
c  wt = 1/(rmsocl^2+1/rmsocl^2)
c  wt = 1/(rmsocl^2+0.1)
c  If rmsocl is very small or very large, then wt is small and that
c  measurement is essentially ignored. This may seem paradoxical,
c  but measurements with extremely small rmsocl are generally found
c  in completely saturated windows in which GFIT has increased CL.

      implicit none
      include "../gfit/ggg_int_params.f"
      include "params.f"
      include "../comn/postproc_params.f"

      integer irow,jtg,j,k,ktg,kt0,lnbc,lloc,idot,lr,
     & iggg,igog,lsf,sflag,
     & npp,ldot,totnit,ntc,mlabel,mval,jval,lunw_nts,
     & lunr_mul,lunr_rl,lunr_col,lunw_rpt,lunw_xsw,lunw_sites,
     & lunw_trp,mcol,ncol,icol,mrow,lnit,lg,ld,tnt,tnd,
     & msite,nsite,jsite,lunw_miss,lunw_vsf,idum,iunder,
     & nauxcol,nlhead,nmiss,naddn,iyrwas,doywas,
c     & nss,
     & nrow,nit,mit,i,iyr,doy,istat,nfound
      parameter (lunw_sites=21)    ! collate_sites.rpt
      parameter (lunr_mul=51)      ! multiggg.sh 
      parameter (lunr_rl=53)       ! runlog
      parameter (lunr_col=55)      ! .col file
      parameter (lunw_xsw=58)      ! output file (.xsw)
      parameter (lunw_rpt=60)      ! .rpt file
      parameter (lunw_trp=61)      ! .trp file (transpose of .xsw)
      parameter (lunw_nts=62)      ! .nts file (negative time step)
      parameter (lunw_miss=65)     ! .missing file
      parameter (lunw_vsf=71)
      parameter (mcol=628)         ! Total number of columns/windows
      parameter (mrow=360000)      ! Max number of output spectra
      parameter (mval=35100000)    ! Max number of values (NROW*NCOL)
      parameter (mlabel=17200)     ! Max # of characters in col labels
      parameter (msite=20000)      ! Max number of observation sites

      character 
     & dl*1,
     & col1*1,
     & apf*2,
     & rlgfile*(mfilepath),
     & version*64,
     & data_fmt_read_rl*256,col_labels_rl*320,
     & ans*1,
     & collabel*(mlabel),
     & sflabel*(mlabel),
     & auxcol*(13*mauxcol),
     & outfile*80,
     & specname_rlg*(nchar),
     & sign(mrow)*1,
     & spectrum(mrow)*(nchar),
     & output_fmt*40,
     & windows(mcol)*12,
     & gggdir(mpath),
     & specname_was*(nchar),
     & addn_lines(maddln)*(mcharhead)

      real*8 asza,
     & zobs,
     & rmin,rmax,vmin,vmax,vemax,vemin,tvsf,tvse,ww,
     & graw,obslat,obslon,opd,
     & r8was,r8year,r8ydiff,year(mrow),
     & sitelat(msite),
     & sitelon(msite),
     & sitealt(msite),
     & wt,twt,trms,lasf,wavtkr,aipl,sia,fvsi,azim,wspd,wdir,osds

      real*4
     & fcen(mcol),width(mcol),
     & airmass,rmsoclpc,fqshift,r2,sg,sfval,
     & vsf,vsf_err,ymiss,small,cl,tilt,cc,zlo,ovcol,
c     & rversion,
     & yaux(mauxcol,mrow),
     & yobs(mval), yerr(mval), qc(mrow),ytot(mrow)
      parameter (ymiss=9.8765e+35,small=1.0E-18)

      real*8 zpdtim,tout,pout,hout,tins,pins,hins,fovi,fovo,amal,
     & snr,zenoff,zoff,zpdwas,max_delta_t,delta_t,zmin

      integer*4 bytepw,ifirst,ilast,possp,ntot(msite),nday(msite),
     & ncount(mrow)

      character yaux_string*600,command_string*256

      logical append_qcflag
      logical append_spectrum_name

      data ntot/msite*0/
      data nday/msite*0/
      data ncount/mrow*0/

      idum=mcolvav  ! Prevent compiler warning (unused paramter)
      idum=mgas     ! Prevent compiler warning (unused paramter)
      idum=mlev     ! Prevent compiler warning (unused paramter)
      idum=mrow_qc  ! Prevent compiler warning (unused paramter)
      idum=mspeci   ! Prevent compiler warning (unused paramter)
      idum=mvmode   ! Prevent compiler warning (unused paramter)
      idum=ncell    ! Prevent compiler warning (unused paramter)

      nsite=0 ! initialize nsite
      sflag=0
      append_qcflag=.false.
      append_spectrum_name=.false.  ! prevent compiler warning
      append_spectrum_name=.true.  ! prevent compiler warning

      call get_ggg_environment(gggdir,dl)
      version=
     &' collate_results          Version 2.08    2020-04-19   GCT,JLL'
      write(6,*) version
      lr=0

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
      elseif(iargc() == 1) then
         call getarg(1, ans)
      else
         stop 'Use: $gggpath/bin/collate_results t/v/l/c/m/n/f/s/r'
      endif

c  Find the number of windows/columns (NCOL)
      open(lunr_mul,file='multiggg.sh',status='old')
      do icol=1,mcol  
121      read(lunr_mul,'(a)',end=99) tabel
         if(tabel(1:1).eq.':') go to 121
         if(lnbc(tabel).le.0) go to 121
      end do  ! icol=1,mcol     !  main loop (over windows)
      read(lunr_mul,*,end=99) tabel
c     write(*,*) 'mcol=',mcol
      stop 'Increase parameter mcol'
 99   close(lunr_mul)
      ncol=icol-1

      open(lunw_vsf,file='collate_results.vsf',status='unknown')
      write(lunw_vsf,*) 2, 6
      write(lunw_vsf,*)' icol  ispec  freq  wid  vsf  vsf_error'
c  Read in the retrieved absorber amounts (YOBS+-YERR)
      open(lunr_mul,file='multiggg.sh',status='old')
      open(lunw_rpt,file='collate_results.rpt',status='unknown')
      open(lunw_nts,file='collate_results.nts',status='unknown')
      write(lunw_rpt,*)2,16
      write(lunw_rpt,*)'iwin  fcen  fcen_error   Nrow   Npp   NIT  %Conv
     &  RMS_Min  RMS_Mean  RMS_Max   VSF_min    VSF_bar    VSF_bar_error
     & VSF_max    &VERR_min    VERR_max'
      do icol=1,ncol     !  main loop (over windows / .col files)
         npp=0
         twt=0.0d0
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
            if(k.eq.6) rlgfile=header_string(:90)    ! GCT 2009-03-04
            if(index(header_string,'runlogs').gt.0)
     &      rlgfile=header_string(:mpath)
         end do
         csformat=header_string(:lnbc(header_string))
c         write(*,'(a)') 'csformat='//csformat

         read(lunr_col,'(a)') command_string
         lsf=index(command_string,'sf=')
         if(lsf.gt.0) then
            sflag=1
            read(command_string(lsf+3:),*) sfval
         else
            sfval=1.0
         endif
         write(sflabel(1+8*(icol-1):),'(f8.3)')sfval

c         write(*,'(a)') 'command_string='//command_string
         read(command_string,*) fcen(icol), width(icol), mit
         read(lunr_col,'(a)')header_string
c         write(*,'(a)') 'header_string='//header_string
         lnit= index(header_string,'Nit')
c         write(*,*) 'lnit= ',lnit

         kt0=index(header_string,' OVC_')
         if(colfile(1:1).eq.'m') then ! Kludge for InSb+InGaAs
           ktg=1+index(header_string,'OVC_'//colfile(2:idot-1))
         else
           ktg=1+index(header_string,' OVC_'//colfile(:idot-1))
         endif
c         write(*,*) 'ktg= ',ktg,' OVC_'//colfile(:idot-1)
         ktg=1+(ktg-kt0)/42
c        ktg=1+(ktg-80)/48
c         write(*,*) 'ktg= ',ktg,' OVC_'//colfile(:idot-1)
c         write(*,*)'colfile,ktg=',colfile(:idot-1),ktg
c         if ( ktg .gt. 1) then
c           call substr(header_string(:ktg-1),cdum,1,nss)
c           ktg=(nss-4)/4
c         endif
         iyrwas=-99999
         doywas=-99999
         zpdwas=-99999.9d0
         irow=0
         rmin=9999999.0
         rmax=0.0
         vmin=+1.E+38
         vmax=-1.E+38
         vemin=+1.E+38
         vemax=-1.E+38
         tvsf=0.0d0
         tvse=0.0d0
         specname_rlg=' '
         specname_was='x'
         nsite=0

         lr=lnbc(rlgfile)
         if(rlgfile(lr-2:lr-2).eq.'o') then
            max_delta_t=0.0004  ! 1.44s (ACE)
         else                 ! Ground-based
c            max_delta_t=0.0025  ! 9.0s 
            max_delta_t=0.0014  ! 5.0s  GCT 2018-08-07
         endif

c  Add spectrum name to output files only on non-MkIV gnd data
         ld=lloc(rlgfile,dl)
         if(rlgfile(lr-2:lr-2).eq.'g' .and.
     &      rlgfile(ld+1:ld+3).ne.'m4_' .and.
     &      rlgfile(ld+1:ld+5).ne.'synth' .and.
     &      rlgfile(ld+1:ld+10).ne.'ll20101005' .and.
     &      rlgfile(ld+1:ld+4).ne.'mkiv' ) then
            append_spectrum_name = .true.
         else
            append_spectrum_name = .false.
         endif
c
c  Read auxilliary measurements from runlog
c         write(*,*) 'runlog=',rlgfile
         open(lunr_rl,file=rlgfile, status='old')   !DG000906
         call read_runlog_header(lunr_rl,data_fmt_read_rl,col_labels_rl)
         r8was=-9999999.9d0
         do while (irow.lt.mrow)  ! Loop over RL recs with different times
            specname_col='='
            read(lunr_col,'(a)',end=24) col_string
            if (lnbc(col_string).le.2) goto 24 ! skip blank line at EOF
c            write(*,'(a)')col_string(:lnbc(col_string))

            read(col_string,csformat)
     &      specname_col(:lnit-3),nit,cl,tilt,cc,fqshift,
     &      sg,zlo,rmsoclpc,zmin,(airmass,ovcol,vsf,vsf_err,jtg=1,ktg)
c            write(*,*)'icol,ktg=',icol,ktg,specname_col(:lnit-3),vsf
            totnit=totnit+nit
            if(nit.lt.mit) ntc=ntc+1  ! Number of Times Converged
            if(rmsoclpc.le.0.0) then
               write(*,*) 'rmsoclpc <= 0  ',colfile,irow,lnit
c               rmsoclpc=0.0001
c               stop 'rmsoclpc <= 0'   ! Commented 2009-03-18
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
            write(lunw_vsf,*)icol,irow+1,fcen(icol),width(icol)/2,
     &      vsf,vsf_err
24          continue

            do while(specname_col.ne.specname_rlg)
               call read_runlog_data_record(lunr_rl,data_fmt_read_rl,
     &         col1,specname_rlg,iyr,doy,zpdtim,obslat,obslon,zobs,
     &         asza,zenoff,azim,osds,opd,fovi,fovo,amal,ifirst,ilast,
     &         graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     &         tout,pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)
c               write(*,*)'read_runlog: irow,istat=',irow+1,istat,
c     &         specname_rlg
c               if(istat.ne.0)write(*,*)'Called readrunlog: istat=',istat
               if(istat.ne.0) go to 14     ! Exit Loop  irow=1,mrow
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
                  if(jsite.eq.9)write(56,*)irow,zobs,iyr,doy,nday(jsite)
               endif

               lg=lnbc(specname_rlg)
               ldot=index(specname_rlg,'.')
               if(ldot.eq.0) ldot=lnbc(specname_rlg)+1
               delta_t=zpdtim-zpdwas
c               write(*,*)lg,ldot,specname_rlg(:lg),delta_t, max_delta_t
c
c  Create speparate entries in the .vsw file for spectra whose
c  ZPD times differ by more than MAX_DELTA_T
c               write(*,*)'iyrwas,iyr=',iyrwas,iyr
c               write(*,*)'doywas,doy=',doywas,doy
               if(iyr.ne.iyrwas .or. doy.ne.doywas .or.
     &         dabs(delta_t).ge.max_delta_t ) then  ! New observation
c
c  The following if-statement shouldn't be necessary. But occasionally
c  you get simultaneous InGaAs/Si scans with very different ZPD times.
c  You don't want them to have separate entries in the .vsw file. So
                  if(specname_rlg(4:lg).ne.specname_was(4:lg)
c                  if(specname_rlg(4:ldot-2).ne.specname_was(4:ldot-2)
c     &           .or. specname_rlg(ldot:).ne.specname_was(ldot:)
     &            ) then
                     irow=irow+1
                     sign(irow)=col1
                     spectrum(irow)=specname_rlg
c  Code works because 2000 was leap year (it won't work for 2100)
                     if(mod(iyr,4).eq.0) then  ! Leap Year
                        r8year=iyr+(doy+zpdtim/24.0d0)/366.0d0
                     else                     ! non-Leap_year
                        r8year=iyr+(doy+zpdtim/24.0d0)/365.0d0
                     endif
                     year(irow)=r8year ! added by DW 20120105
                     write(yaux_string,*)doy+zpdtim/24.d0,zpdtim,irow,
     &               obslat,obslon,zobs,zmin,asza+zenoff,azim,osds,opd,
     &               fovi,amal,graw,tins,pins,tout,pout,
     &               hout,sia,fvsi,wspd,wdir
                     read(yaux_string,*) (yaux(j,irow),j=2,24)

c                     yaux(2,irow)=doy+zpdtim/24.0d0
c                     yaux(3,irow)=zpdtim
c                     yaux(4,irow)=irow
c                     yaux(5,irow)=obslat
c                     yaux(6,irow)=obslon
c                     yaux(7,irow)=zobs
c                     yaux(8,irow)=zmin
c                     yaux(9,irow)=asza+zenoff
c                     yaux(10,irow)=azim
c                     yaux(11,irow)=osds
c                     yaux(12,irow)=opd
c                     yaux(13,irow)=sqrt(fovi**2+amal**2)
c                     yaux(14,irow)=graw
c                     yaux(15,irow)=tins
c                     yaux(16,irow)=pins
c                     yaux(17,irow)=tout
c                     yaux(18,irow)=pout
c                     yaux(19,irow)=hout
c                     yaux(20,irow)=sia
c                     yaux(21,irow)=fvsi
c                     yaux(22,irow)=wspd
c                     yaux(23,irow)=wdir
                     jval=jval+ncol   !  jval=icol+ncol*(irow-1)
                     yobs(jval)=ymiss
                     yerr(jval)=ymiss
                  endif  !   if(specname_rlg.ne.specname_was
               else
c                  write(*,*)'not new obs '//specname_rlg(:20),delta_t
               endif   ! if(iyr.ne.iyrwas .or. doy.ne.doywas 
               iyrwas=iyr
               doywas=doy
               zpdwas=zpdtim
               specname_was=specname_rlg
c               write(*,*)'irow,specname_was: ',irow,specname_was
            end do     ! while(specname_col.ne.specname_rlg)
        
c            irow=irow+1
            sign(irow)=col1
            spectrum(irow)=specname_rlg
            ncount(irow)=ncount(irow)+1
c            write(*,*) icol,irow,ncount(irow),spectrum(irow)
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
            year(irow)=r8year
c            yaux(1,irow)=r8year

            write(yaux_string,*)doy+zpdtim/24.d0,zpdtim,irow,
     &      obslat,obslon,zobs,zmin,asza+zenoff,azim,osds,opd,
     &      fovi,amal,graw,tins,pins,tout,pout,
     &      hout,sia,fvsi,wspd,wdir
            read(yaux_string,*) (yaux(j,irow),j=2,23)
c               write(*,*)'2:asza,zenoff=',asza,zenoff,asza+zenoff
c               write(*,*)'yaux(9)=',yaux(9,irow)

c            yaux(2,irow)=doy+zpdtim/24.0d0
c            yaux(3,irow)=zpdtim
c            yaux(4,irow)=irow
c            yaux(5,irow)=obslat
c            yaux(6,irow)=obslon
c            yaux(7,irow)=zobs
c            yaux(8,irow)=zmin
c            yaux(9,irow)=asza+zenoff
c            yaux(10,irow)=azim
c            yaux(11,irow)=osds
c            yaux(12,irow)=opd
c            yaux(13,irow)=sqrt(fovi**2+amal**2)
c            yaux(14,irow)=graw
c            yaux(15,irow)=tins
c            yaux(16,irow)=pins
c            yaux(17,irow)=tout
c            yaux(18,irow)=pout
c            yaux(19,irow)=hout
c            yaux(20,irow)=sia
c            yaux(21,irow)=fvsi
c            yaux(22,irow)=wspd
c            yaux(23,irow)=wdir
c
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
               yobs(jval)=1000*fqshift/fcen(icol)
               yerr(jval)=vsf_err
c              yerr(jval)=vsf_err*abs(yobs(icol,irow))
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
            r2=3.0*abs(rmsoclpc)  ! rmsoclps typically ~ 0.33%
c            wt= 1./(r2+1.0/r2)
            wt= 1./(r2+0.1)
            twt=twt+wt
            trms=trms+wt*rmsoclpc
         end do   !  irow.lt.mrow  !  Loop over spectra
         stop 'Increase parameter MROW'
14       nrow=irow ! the number of records in the .col and .vsw files
         if(ncol*nrow.gt.mval) then
            write(*,*)ncol,nrow
            write(*,*)'Increase parameter MVAL to ',nrow*ncol
            stop 'collate'
         endif
         close(lunr_rl)
         close(lunr_col)
         write(lunw_rpt,'(i3,f9.2,f8.2,2i8,2f7.2,3f9.4,6(1pe12.4))')
     &   icol,fcen(icol),width(icol)/2,
     &   npp,nrow,float(totnit)/npp,100*float(ntc)/npp,
     &   rmin,trms/twt,rmax,vmin,tvsf/tvse,sqrt(tvse*nrow)/tvsf,
     &   vmax,vemin,vemax
      end do  ! icol=1,mcol     !  main loop (over windows)
      close(lunr_mul)
      close(lunw_rpt)
      close(lunw_nts)
c====================================================================
c      do k=lr,1,-1
c         if(ichar(rlgfile(k:k)) .eq. 92) go to 101  ! backslash
c         if(ichar(rlgfile(k:k)) .eq. 47) go to 101  ! forward slash
c      end do
c101   outfile=rlgfile(k+1:lr-3)//ans//'sw'
      outfile=rlgfile(lloc(rlgfile,dl)+1:lr-3)//ans//'sw'
c==================================================================
      auxcol=
     &'      year         day          hour         run    '//
     &'      lat          long         zobs         zmin   '//
     &'      solzen       azim         osds         opd    '//
     &'      fovi         amal         graw         tins   '//
     &'      pins         tout         pout         hout   '//
     &'      sia          fvsi         wspd         wdir   '
c====================================================================
      call substr(auxcol,cdum,1,nauxcol)
      if (append_qcflag) then
         auxcol=auxcol(:lnbc(auxcol))//'      qcflag '
         nauxcol=nauxcol+1
         call generate_qc_flag(nrow,spectrum,qc)
         do irow=1, nrow
            yaux(nauxcol,irow)=qc(irow)
         enddo
      endif
c====================================================================
c      call substr(auxcol,cdum,1,nauxcol)
c      if(nauxcol.gt.mauxcol) then
c          write(*,*)' mauxcol, nauxcol = ',mauxcol,nauxcol
c          stop 'increase parameter mauxcol'
c      endif

      output_fmt='(a1,f13.8,NNf13.5,MMMM(1pe12.4))'
      write(output_fmt(11:12),'(i2.2)') nauxcol-1
      write(output_fmt(19:22),'(i4.4)') 2*ncol
      if (append_spectrum_name) then
         auxcol='  spectrum   '//auxcol(:13*nauxcol)
         output_fmt='(a57,'//output_fmt(2:)
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
      if(sflag .eq. 1) then
        addn_lines(4) = 'sf='//sflabel(:lnbc(sflabel))
        naddn = 4
      endif
      call write_postproc_header(lunw_xsw, nauxcol+2*ncol, nrow, 
     & nauxcol, dble(ymiss), output_fmt, addn_lines, naddn, 0)
      write(lunw_xsw,'(a)') auxcol(:lnbc(auxcol))//'  '//
     &   collabel(:lnbc(collabel)+1)

c  Note that the SIGN array and the following IF statement are merely
c  to support both the new and the old runlog formats.
      open(lunw_miss,file='collate_results.missing',status='unknown')
      nfound=0
      nmiss=0
      jval=0
      do irow=1,nrow
         do k=1,ncol
            jval=jval+1     !    jval=k+ncol*(irow-1)
c            if( yobs(jval).eq.ymiss .and.
c     &          yerr(jval).eq.ymiss) then
            if( abs(yobs(jval)-ymiss).lt.small*ymiss .and.
     &          abs(yerr(jval)-ymiss).lt.small*ymiss) then
               nmiss=nmiss+1
               write(lunw_miss,*)'Missing: ',windows(k),
     &       '  '//spectrum(irow)
            else
               nfound=nfound+1
            endif
         end do
         if (append_spectrum_name) then
            write(lunw_xsw,output_fmt) spectrum(irow),sign(irow),
     &      year(irow),(yaux(k,irow),k=2,nauxcol-1),
     &      (yobs(k+ncol*(irow-1)),yerr(k+ncol*(irow-1)),k=1,ncol)
         else
            write(lunw_xsw,output_fmt) sign(irow),year(irow),
     &      (yaux(k,irow),k=2,nauxcol),
     &      (yobs(k+ncol*(irow-1)),yerr(k+ncol*(irow-1)),k=1,ncol)
         endif

      end do  !  irow=1,nrow
      close(lunw_miss)
      close(lunw_xsw)
c====================================================================
      do j=1,nrow
         ytot(j)=0.0
      end do
c  Write another output file containing the transpose of data matrix
      open(lunw_trp,file=outfile(:lnbc(outfile))//'.transpose',
     & status='unknown')
      write(lunw_trp,*)4,nrow+3
      write(lunw_trp,*)'Missing:',ymiss
      write(lunw_trp,'(a24,600000a36)')' Window fcen fcen_error ',
     & (spectrum(j),j=1,nrow)
      write(lunw_trp,'(a24,600000i4)')' Window fcen fcen_error ',
     & (j,j=1,nrow)
      do k=1,ncol
         write(lunw_trp,*)k,fcen(k),1+width(k)/2,
     &  (yobs(k+ncol*(j-1)),j=1,nrow)
         do j=1,nrow
            ytot(j)=ytot(j)+yobs(k+ncol*(j-1))
         end do
      end do
      close(lunw_trp)

c     write(67,*) 2,3
c     write(67,*) ' j zmin rms'
c     do j=1,nrow
c        write(67,*) j, yaux(8,j), ytot(j)/ncol
c     end do
c====================================================================
      write(6,*)
      write(6,*) outfile(:lnbc(outfile))//' contains:'
      write(6,'(i8,a)') nauxcol,' auxiliary columns'
      write(6,'(i8,a)') ncol,' data column pairs (value + error)'
      write(6,'(i8,a)') 2*ncol+nauxcol,' total columns'
      write(6,'(i8,a)') nrow,' data rows'
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
      write(lunw_sites,*)' #   Nobs   Nday   Lat   Long    Alt'
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
