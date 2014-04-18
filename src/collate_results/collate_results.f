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
c     runlog.xsw   Spreadsheet of individual window values for each spectrum

      implicit none
      include "../ggg_int_params.f"
      include "params.f"

      integer i1,idot,irow,jtg,j,k,ktg,lnbc,fnbc,fbc,lloc,
     & npp,ldot,l2,l3,l4,totnit,ntc,mlabel,mval,jval,lunw_nts,
     & lr,lunr_mul,lunr_rl,lunr_col,lunw_rpt,lunw_xsw,luns,
     & mcol,ncol,icol,mrow,lnit,lg,ld,
     & msite,nsite,jsite,lunw_miss,
     & nauxcol,nlhead,nmiss,iyrwas,doywas,
     & nrow,nss,nit,mit,i,iyr,doy,istat,nfound
      parameter (luns=21)           ! collate_sites.rpt
      parameter (lunr_mul=51)       ! multiggg.sh 
      parameter (lunr_rl=53)        ! runlog
      parameter (lunr_col=55)       ! .col file
      parameter (lunw_xsw=58)       ! output file (.xsw)
      parameter (lunw_rpt=60)       ! .rpt file
      parameter (lunw_nts=62)       ! .nts file (negative time step)
      parameter (lunw_miss=65)      ! .nts file (negative time step)
      parameter (mcol=600)          ! Total number of columns/windows
      parameter (mrow=360000)       ! Max number of output records/spectra
      parameter (mval=35100000)     ! Max number of values (NROW * NCOL)
      parameter (mlabel=16000)      ! Max # of characters in column labels
      parameter (msite=1500)          ! Max number of observation sites

      character 
     & dl*1,
     & col1*1,
     & apf*2,
     & rlgfile*(mfilepath),
     & version*64,
     & data_fmt_read_rl*256,col_labels_rl*320,
     & ans*1,
     & collabel*(mlabel),
     & auxcol*(13*mauxcol),
     & outfile*80,
     & specname_grl*(nchar),
     & sign(mrow)*1,
     & spectrum(mrow)*(nchar),
     & output_fmt*40,
     & windows(mcol)*10,
     & specname_gwas*(nchar)

      real*8 airmass,asza,cl,tilt,cc,zlo,
     & fcen(mcol),width(mcol),
     & zobs,
     & rmin,rmax,vmin,vmax,vemax,vemin,tvsf,tvse,ww,
     & fqshift,graw,obslat,obslon,opd,ovcol,rmsfit,
     & r8was,r8year,r8ydiff,year(mrow),
     & sitelat(msite),
     & sitelon(msite),
     & sitealt(msite),
     & trms,lasf,wavtkr,aipl,sia,fvsi,azim,wspd,wdir,osds

      real*4
     & vsf,vsf_err,ymiss,
c     & rversion,
     & yaux(mauxcol,mrow),
     & yobs(mval), yerr(mval), qc(mrow),ytot(mrow)
      parameter (ymiss=9.8765e+35)

      real*8 zpdtim,tout,pout,hout,tins,pins,hins,fovi,fovo,amal,
     & snr,zenoff,zoff,sg,zpdwas,max_delta_t,delta_t,zmin

      integer bytepw,ifirst,ilast,possp,ntot(msite)

      logical append_qcflag
      logical append_spectrum_name

      data ntot/msite*0/

      nsite=0 ! initialize nsite
      append_qcflag=.false.
      append_spectrum_name=.true.  ! prevent compiler warning
      append_spectrum_name=.false.  ! prevent compiler warning

      version=
     &' collate_results             Version 1.71     2013-10-18   GCT'
      write(6,*) version
      lr=0

c  Initialize character arrays (Necessary for the G77 compiler).
      do i=1,mlabel
         collabel(i:i)=' '
      end do

      if (iargc() == 0) then
         write(6,'(a)')
     $' vsf [t], vertical column [v], los column [l],',
     $' continuum [c], tilt [m], cc [n],  frequency shift [f],', 
     &'  solar-gas frequency shift [s], rms [r]: '
         read(5,*) ans
      elseif(iargc() == 1) then
         call getarg(1, ans)
      else
         stop 'Use: $gggpath/bin/collate_results t (or v/l/c/m/n/f/s/r)'
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

c  Read in the retrieved absorber amounts (YOBS+-YERR)
      open(lunr_mul,file='multiggg.sh',status='old')
      open(lunw_rpt,file='collate_results.rpt',status='unknown')
      open(lunw_nts,file='collate_results.nts',status='unknown')
      write(lunw_rpt,*)2,15
      write(lunw_rpt,*)' iwin     fcen  fcen_error   Nrow     Npp
     &   NIT  %Conv   RMS_Min  RMS_Mean  RMS_Max   VSF_min     VSF_bar  
     &   VSF_max    VERR_min   VERR_max'
      do icol=1,ncol     !  main loop (over windows)
        npp=0
        trms=0.0d0
        totnit=0
        ntc=0
        jval=icol-ncol
135     read(lunr_mul,'(a)') tabel
        if(tabel(1:1).eq.':') go to 135
         if(lnbc(tabel).le.0) go to 135
        colfile=tabel(index(tabel,' ')+1:index(tabel,'.ggg'))//'col'
        write(*,*) 'Opening '//colfile(:64)
        open(lunr_col,file=colfile,status='old') ! .col file
        idot=index(colfile,'.')
        i1=index(colfile(:idot),'_')
        if(i1.eq.0) i1=index(colfile(:idot),'^')   !  new format
        if(i1.eq.0) i1=index(colfile(:idot),'-')   !  old format
        collabel=collabel(:lnbc(collabel)+2)//colfile(:idot-1)//' '//
     $  colfile(:idot-1)//'_error'
        windows(icol)=colfile(:idot-1)
c
c  Read header lines of .col file and locate column containing "OVC_gas".
c  in order to read data from appropriate target gas.
        read(lunr_col,*) nlhead
        read(lunr_col,'(a)') gfit_version
        read(lunr_col,'(a)') gsetup_version
        do k=4,nlhead-2
           if(nlhead.ge.23) read(lunr_col,'(34x,a)')header_string
           if(nlhead.eq.22) read(lunr_col,'(34x,a)')header_string
           if(nlhead.eq.20) read(lunr_col,'(34x,a)')header_string
           if(nlhead.eq.21) read(lunr_col,'(a)')header_string
           if(k.eq.6) rlgfile=header_string(:90)    ! GCT 2009-03-04
           if(index(header_string,'runlogs').gt.0)
     &     rlgfile=header_string(:mpath)
        end do
        dl=rlgfile(1:1)
        csformat=header_string(:lnbc(header_string))
c        write(*,'(a)') 'csformat='//csformat

        read(lunr_col,'(a)') header_string
c        write(*,'(a)') 'header_string='//header_string
        read(header_string,*) fcen(icol), width(icol), mit
c        if(icol.gt.1) write(44,*)fcen(icol-1),fcen(icol),
c     &  fcen(icol-1)+width(icol-1)/2,fcen(icol)-width(icol)/2,
c     &  fcen(icol)-width(icol)/2-fcen(icol-1)-width(icol-1)/2
        read(lunr_col,'(a)')header_string
        lnit= index(header_string,'Nit')

c        write(*,*) header_string
c        write(*,*) colfile(:i1-1)
        ktg=1+index(header_string,' OVC_'//colfile(:i1-1))
c        write(*,*)'ktg=',ktg
        if ( ktg .gt. 1) then
          call substr(header_string(:ktg-1),cdum,1,nss)
          ktg=(nss-4)/4
        endif
c        write(*,*)'ktg=',ktg
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
        specname_grl=' '
        specname_gwas='x'
        nsite=0

        lr=lnbc(rlgfile)
        if(rlgfile(lr-2:lr-2).eq.'o') then
           max_delta_t=0.0004  ! 1.44s (ACE)
        else
           max_delta_t=0.0025  ! 9.0s (ground-based TCCON)
        endif

c  Add spectrum name to output files only on non-MkIV gnd data
        ld=lloc(rlgfile,dl)
        if(rlgfile(lr-2:lr-2).eq.'g' .and.
     &     rlgfile(ld+1:ld+3).ne.'m4_' .and.
     &     rlgfile(ld+1:ld+10).ne.'ll20101005' .and.
     &     rlgfile(ld+1:ld+4).ne.'mkiv' ) then
           append_spectrum_name = .true.
        else
           append_spectrum_name = .false.
        endif
c
c  Read auxilliary measurements from runlog
c        write(*,*) 'runlog=',rlgfile
        open(lunr_rl,file=rlgfile, status='old')   !DG000906
        call read_runlog_header(lunr_rl,data_fmt_read_rl,col_labels_rl)
c        read(lunr_rl,*) nlhead,nn
c        do i=2,nlhead
c           read(lunr_rl,*)
c        end do
        r8was=-9999999.9d0
        do while (irow.lt.mrow)   !  Loop over runlog records with different times
           specname_col='='
           read(lunr_col,'(a)',end=24) col_string
           if ( lnbc(col_string) .le. 2 ) go to 24  ! skip blank line at EOF
c           write(*,'(a)')col_string(:lnbc(col_string))
           l2=fbc(col_string(2:))+1       ! First space following spectrum name
           l3=fnbc(col_string(l2:))+l2-1  ! First character of NIT
           l4=fbc(col_string(l3:))+l3-1   ! First space following NIT
c           write(*,*)l2,l3,l4

c           if (index(gfit_version,'2.40.2') .ne. 0) then !  old col file format 
c               csformat='(1x,a21,i2,1x,f5.3,3(1x,f4.1),1x,f5.3,'
c     &         //'1x,f6.4,f7.3,1x,9(0pf7.3,1pe10.3,0pf9.4,1pe8.1))'
c           elseif (index(gfit_version,'4.8.') .ne. 0) then !  new .col file format
c               write(csformat,'(a,i2.2,a)')'(1x,a',l4-5,
c     &    ',i3,f6.3,4f5.1,f6.3,f7.4,f8.3,15(f7.3,e11.4,f9.4,e8.1))'
c           else                                  ! assume previous .col file format 
c               write(csformat,'(a,i2.2,a)')'(1x,a',l4-5,
c     &     ',i3,f6.3,3f5.1,f6.3,f7.4,f8.3,15(f7.3,e11.4,f9.4,e8.1))'
c           endif

c            lv=index(gfit_version, 'Version ')
c            read(gfit_version(lv+7:),*) rversion
c            if (rversion .ge. 4.8) then
           if (index(gfit_version,'4.8') .gt. 0 .or.
     &       index(gfit_version,'4.9') .gt. 0 .or.
     &       index(gfit_version,'4.10') .gt. 0 .or.
     &       index(gfit_version,'4.11') .gt. 0 .or.
     &       index(gfit_version,'4.12') .gt. 0 .or.
     &       index(gfit_version,'4.13') .gt. 0 .or.
     &       index(gfit_version,'4.14') .gt. 0 .or.
     &       index(gfit_version,'4.20') .gt. 0 .or.
     &       index(gfit_version,'4.21') .gt. 0 .or.
     &       index(gfit_version,'4.22') .gt. 0 .or.
     &       index(gfit_version,'4.23') .gt. 0 .or.
     &       index(gfit_version,'4.24') .gt. 0 .or.
     &       index(gfit_version,'4.25') .gt. 0 .or.
     &       index(gfit_version,'4.26') .gt. 0 .or.
     &       index(gfit_version,'4.27') .gt. 0 .or.
     &       index(gfit_version,'4.28') .gt. 0 .or.
     &       index(gfit_version,'4.29') .gt. 0 .or.
     &       index(gfit_version,'4.30') .gt. 0 .or.
     &       index(gfit_version,'4.31') .gt. 0 .or.
     &       index(gfit_version,'4.32') .gt. 0 .or.
     &       index(gfit_version,'4.33') .gt. 0 .or.
     &       index(gfit_version,'4.34') .gt. 0 .or.
     &       index(gfit_version,'4.35') .gt. 0 .or.
     &       index(gfit_version,'4.36') .gt. 0 .or.
     &       index(gfit_version,'4.37') .gt. 0 
     &       )    then   !  .col file includes CC
             read(col_string,csformat)
     &       specname_col(:lnit-3),nit,cl,tilt,cc,fqshift,
     &       sg,zlo,rmsfit,zmin,(airmass,ovcol,vsf,vsf_err,jtg=1,ktg)
           else                                          !  no CC
             read(col_string,csformat)
     &       specname_col(:lnit-3),nit,cl,tilt,fqshift,
     &       sg,zlo,rmsfit,zmin,(airmass,ovcol,vsf,vsf_err,jtg=1,ktg)
           endif
c           write(*,*) specname_col,nit,cl,tilt,cc,fqshift,
c     &     sg,zlo,rmsfit,zmin,airmass,ovcol,vsf,vsf_err
           totnit=totnit+nit
           if(nit.lt.mit) ntc=ntc+1  ! Number of Times Converged
           if(rmsfit.le.0.0) then
             write(*,*) 'rmsfit <= 0  ',colfile,irow
             write(*,*)col_string(:72)
             rmsfit=0.0001
c             stop 'rmsfit <= 0'   ! Commented 2009-03-18
           endif
           if(rmsfit.gt.rmax) rmax=rmsfit
           if(rmsfit.lt.rmin) rmin=rmsfit
           if(vsf.lt.vmin) vmin=vsf
           if(vsf.gt.vmax) vmax=vsf
           if(vsf_err.lt.vemin) vemin=vsf_err
           if(vsf_err.gt.vemax) vemax=vsf_err
           ww=1.0d0/vsf_err**2
           tvse=tvse+ww
           tvsf=tvsf+vsf*ww
24         continue

           do while(specname_col.ne.specname_grl)
              call read_runlog_data_record(lunr_rl,data_fmt_read_rl,
     &        col1,specname_grl,iyr,doy,zpdtim,obslat,obslon,zobs,
     &        asza,zenoff,azim,osds,opd,fovi,fovo,amal,ifirst,ilast,
     &        graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     &        tout,pout,hout,sia,fvsi,wspd,wdir,lasf,wavtkr,aipl,istat)
c              if(istat.ne.0) write(*,*) 'istat=',istat
c              write(*,*) '...called read_runlog.'
              if(istat.ne.0) go to 14     ! Exit Loop  irow=1,mrow
              if(col1.eq.':' .or. col1.eq.';') cycle
c              if(wspd.eq.-1.0) wspd=ymiss
c              if(wdir.eq.-1.0) wdir=ymiss

              do jsite=1,nsite
                 if(obslat.eq.sitelat(jsite)
     &           .and.obslon.eq.sitelon(jsite)
     &           .and. zobs.eq.sitealt(jsite)) go to 66
              end do
              nsite=nsite+1
              sitelat(jsite)=obslat
              sitelon(jsite)=obslon
              sitealt(jsite)=zobs
66            ntot(jsite)=ntot(jsite)+1

              lg=lnbc(specname_grl)
              ldot=index(specname_grl,'.')
              delta_t=zpdtim-zpdwas
c              write(*,*)specname_grl(:lg),delta_t, max_delta_t
c
c  Create speparate entries in the .vsw file for spectra whose
c  ZPD times differ by more than MAX_DELTA_T
c              write(*,*)'iyr=',iyrwas,iyr
c              write(*,*)'doy=',doywas,doy
              if(iyr.ne.iyrwas .or. doy.ne.doywas .or.
     &        dabs(delta_t).ge.max_delta_t ) then  ! New observation
c
c  The following if-statement shouldn't be necessary. But occasionally
c  you get simultaneous InGaAs/Si scans with very different ZPD times.
c  You don't want them to have separate entries in the .vsw file. So....
c                 write(*,*)specname_grl,specname_gwas
                 if(specname_grl(4:ldot-2).ne.specname_gwas(4:ldot-2)
     &           .or.  specname_grl(ldot:).ne.specname_gwas(ldot:)) then

                    irow=irow+1
                    sign(irow)=col1
                    spectrum(irow)=specname_grl
                    r8year=iyr+(doy+zpdtim/24.0d0)/366.0d0
                    year(irow)=r8year ! added by DW 20120105
                    yaux(2,irow)=doy+zpdtim/24.0d0
                    yaux(3,irow)=zpdtim
                    yaux(4,irow)=irow
                    yaux(5,irow)=obslat
                    yaux(6,irow)=obslon
                    yaux(7,irow)=zobs
                    yaux(8,irow)=zmin
                    yaux(9,irow)=asza+zenoff
                    yaux(10,irow)=azim
                    yaux(11,irow)=osds
                    yaux(12,irow)=opd
                    yaux(13,irow)=sqrt(fovi**2+amal**2)
                    yaux(14,irow)=graw
                    yaux(15,irow)=tins
                    yaux(16,irow)=pins
                    yaux(17,irow)=tout
                    yaux(18,irow)=pout
                    yaux(19,irow)=hout
                    yaux(20,irow)=sia
                    yaux(21,irow)=fvsi
                    yaux(22,irow)=wspd
                    yaux(23,irow)=wdir
                    jval=jval+ncol   !  jval=icol+ncol*(irow-1)
                    yobs(jval)=ymiss
                    yerr(jval)=ymiss
                  endif
              endif
              iyrwas=iyr
              doywas=doy
              zpdwas=zpdtim
              specname_gwas=specname_grl
           end do     ! while(specname_col.ne.specname_grl)
        
           sign(irow)=col1
           spectrum(irow)=specname_grl
           r8year=iyr+(doy+zpdtim/24.0d0)/366.0d0
           r8ydiff=r8year-r8was
c  Report negative time-steps in the runlogs times.
c  Need only do this for the first window (avoid repetition).
           if( r8ydiff .lt. -0.00000001d0 .and. icol.eq.1)
     &     write(lunw_nts,'(a,a,2f12.6)')
     &  '  Negative time step (runlog unsorted?) ',
     *     specname_grl,r8was,r8year
           r8was=r8year
           year(irow)=r8year
c           yaux(1,irow)=r8year
           yaux(2,irow)=doy+zpdtim/24.0d0
           yaux(3,irow)=zpdtim
           yaux(4,irow)=irow
           yaux(5,irow)=obslat
           yaux(6,irow)=obslon
           yaux(7,irow)=zobs
           yaux(8,irow)=zmin
           yaux(9,irow)=asza+zenoff
           yaux(10,irow)=azim
           yaux(11,irow)=osds
           yaux(12,irow)=opd
           yaux(13,irow)=sqrt(fovi**2+amal**2)
           yaux(14,irow)=graw
           yaux(15,irow)=tins
           yaux(16,irow)=pins
           yaux(17,irow)=tout
           yaux(18,irow)=pout
           yaux(19,irow)=hout
           yaux(20,irow)=sia
           yaux(21,irow)=fvsi
           yaux(22,irow)=wspd
           yaux(23,irow)=wdir
c
           if(ans.eq.'t') then
               yobs(jval)=vsf
               yerr(jval)=vsf_err
           elseif(ans.eq.'v') then
               yobs(jval)=vsf*ovcol
               yerr(jval)=vsf_err*ovcol
           elseif(ans.eq.'l') then
               yobs(jval)=vsf*ovcol*airmass
               yerr(jval)=vsf_err*ovcol*airmass
               if(yerr(jval).gt.3.4028E+38) yerr(jval)=3.4028E+38
               if(yerr(jval).lt.1.175E-38) yerr(jval)=1.175E-38
           elseif(ans.eq.'f') then
               yobs(jval)=1000*fqshift/fcen(icol)
               yerr(jval)=vsf_err
c             yerr(jval)=vsf_err*abs(yobs(icol,irow))
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
               yobs(jval)=rmsfit
               yerr(jval)=vsf_err*rmsfit
           elseif(ans.eq.'c') then
               yobs(jval)=cl
               yerr(jval)=cl*rmsfit
           else
               stop 'unknown option'
           endif
           npp=npp+1
           trms=trms+1.0d0/rmsfit**2
c           write(*,*)irow,npp,specname_col
        end do   !  irow.lt.mrow  !  Loop over spectra
        stop 'Increase parameter MROW'
14      nrow=irow ! the number of records in the .col and .vsw files
        if(ncol*nrow.gt.mval) then
           write(*,*)ncol,nrow
           write(*,*)'Increase parameter MVAL to ',nrow*ncol
           stop 'collate'
        endif
        close(lunr_rl)
        close(lunr_col)
        write(lunw_rpt,'(i3,2f8.2,2i8,2f7.2,3f9.4,5e12.3)')icol,
     & fcen(icol),width(icol)/2,
     & npp,nrow,float(totnit)/npp,100*float(ntc)/npp,
     & rmin,1/sqrt(trms/npp),rmax,vmin,tvsf/tvse,vmax,vemin,vemax
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
     &'      asza         azim         osds         opd    '//
     &'      fovi         graw         tins         pins   '//
     &'      tout         pout         hout         sia    '//
     &'      fvsi         wspd         wdir   '
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
         auxcol='  spectrum   '//auxcol(:12*mauxcol)
         output_fmt='(a57,'//output_fmt(2:)
         nauxcol=nauxcol+1
      endif
      write(*,*) ' Output format = '//output_fmt

c  Write out all analyzed abundances to the .?sw disk file.
      open(lunw_xsw,file=outfile,status='unknown')
      write(lunw_xsw,'(i2,i4,i7,i4)')7,nauxcol+2*ncol,nrow,nauxcol
      write(lunw_xsw,'(a)') version(:lnbc(version))
      write(lunw_xsw,'(a)') gfit_version(:lnbc(gfit_version))
      write(lunw_xsw,'(a)') gsetup_version(:lnbc(gsetup_version))
      write(lunw_xsw,'(a8,1pe12.4)')'missing:', ymiss
      write(lunw_xsw,'(a)') 'format='//output_fmt
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
             if( yobs(jval).eq.ymiss .and. yerr(jval).eq.ymiss ) then
          write(lunw_miss,*)'Missing: ', windows(k),'  '//spectrum(irow)
                nmiss=nmiss+1
             else
                nfound=nfound+1
             endif
          end do
          if (append_spectrum_name) then
             write(lunw_xsw,output_fmt) spectrum(irow),sign(irow),
     &       year(irow),(yaux(k,irow),k=2,nauxcol-1),
     &       (yobs(k+ncol*(irow-1)),yerr(k+ncol*(irow-1)),k=1,ncol)
          else
             write(lunw_xsw,output_fmt) sign(irow),year(irow),
     &       (yaux(k,irow),k=2,nauxcol),
     &       (yobs(k+ncol*(irow-1)),yerr(k+ncol*(irow-1)),k=1,ncol)
          endif

      end do  !  irow=1,nrow
      close(lunw_miss)
      close(lunw_xsw)
c====================================================================
      do j=1,nrow
         ytot(j)=0.0
      end do
c  Write another output file containing the transpose of the data matrix
c     write(66,*)2,nrow+3
c     write(66,'(a24,20a20)')' Window fcen fcen_error ',
c    & (spectrum(j),j=1,nrow)
      do k=1,ncol
c        write(66,*)k,fcen(k),width(k)/2,(yobs(k+ncol*(j-1)),j=1,nrow)
         do j=1,nrow
            ytot(j)=ytot(j)+yobs(k+ncol*(j-1))
         end do
      end do
c     write(67,*) 2,3
c     write(67,*) ' j zmin rms'
c     do j=1,nrow
c        write(67,*) j, yaux(8,j), ytot(j)/ncol
c     end do
c====================================================================
      write(6,*)
      write(6,*) outfile(:lnbc(outfile))//' contains:'
      write(6,'(i8,a)') nauxcol,' auxiliary columns'
      write(6,'(i8,a)') ncol,' data parameters (each value + error)'
      write(6,'(i8,a)') 2*ncol+nauxcol,' total columns'
      write(6,'(i8,a)') nrow,' data rows'
      write(6,'(i8,a,f5.1,a1)') nfound,' found values   = ',
     & 100*float(nfound)/(nfound+nmiss),'%'
      write(6,'(i8,a,f5.1,a30)') nmiss,' missing values = ',
     & 100*float(nmiss)/(nfound+nmiss),'%  See collate_results.missing'
c====================================================================
c  Report site information
      open(luns,file='collate_sites.rpt',status='unknown')
      write(luns,*)2,5
      write(luns,*)' #  Nobs  Lat  Long  Alt'
      do jsite=1,nsite
         write(luns,'(i3,i8,3f12.3)') jsite,ntot(jsite)/ncol,
     &   sitelat(jsite),sitelon(jsite),sitealt(jsite)
      end do
      close(luns)
      stop
      end
