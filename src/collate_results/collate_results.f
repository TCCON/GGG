c  Program COLLATE
c
c  Reads the .col output files (gas_1234.runlog.col) produced by GFIT,
c  collates them with auxiliary data from the runlog, and then writes
c  the information to a spreadsheet-format XYPLOT-readable output file
c  named  runlog.tsw, runlog.csw, etc.
c
c  INPUT FILES:
c     multiggg.bat         batch file containing names of .col files
c     runlog.xrl           runlog file
c     gas_1234.runlog.col  files containing column amounts
c
c OUTPUT FILES:
c     runlog.xsw   Spreadsheet of individual window values for each spectrum

      implicit none
      integer i1,idot,irow,j,jtg,k,ktg,lnbc,fnbc,fbc,
     & npp,lrl,l2,l3,l4,totnit,ntc,mlabel,mval,jval,
     &lr,lunm,lunr,lunc,lunv,lunw,mcol,ncol,icol,
     &mrow,nauxcol,nlhead,nmiss,iyrwas,doywas,
     &nrow,nss,nit,mit,i,iyr,doy,istat,nfound
      parameter (lunm=11)      ! multiggg.bat
      parameter (lunr=12)      ! runlog
      parameter (lunc=14)      ! .col file
      parameter (lunw=15)      ! output file (.xsw)
      parameter (lunv=16)      ! .rpt file
      parameter (mcol=600)     ! Total number of columns/windows
      parameter (mrow=160000)  ! Max number of output records/spectra
      parameter (mval=3200000) ! Max number of values (NROW * NCOL)
      parameter (nauxcol=19)   ! Number of auxiliary parameters/columns
      parameter (mlabel=16000) ! Max Number of column lable characters

      character ans*1,apf*2,cdum*20,colabel*500,
     & gfit_version*80,gsetup_version*80,col_string*500,
     & csformat*80,collabel*(mlabel),outfile*20,col1*1,
     & spname_rl*21,specname_col*21,runlog*72,sign(mrow)*1,
     & spectrum(mrow)*21,tabel*80,
     & colfile*40,collate_version*48,window(mcol)*10,spname_rlwas*21

      real*8 airmass,asza,cl,tilt,zlo,fcen,width,zobs,rmin,rmax,
     & fqshift,graw,obslat,obslon,opd,ovcol,rmsfit,
     & r8was,r8year,r8ydiff,year(mrow),
     & trms,lasf,wavtkr,sia,sis,aipl

      real*4
     & vsf,vsf_err,ymiss,
     & yaux(nauxcol,mrow),
     & yobs(mval), yerr(mval)
      parameter (ymiss=9.9999e+29)

      real*8 zpdtim,tout,pout,hout,tins,pins,hins,fovi,fovo,amal,
     & snr,zenoff,zoff,sg,zpdwas,max_delta_t,delta_t,zmin
      integer bytepw,ifirst,ilast,possp

      collate_version=' collate_results  Version 1.0.0  2008-10-03  GCT'
      write(6,*) collate_version

c  Initialize character arrays (Necessary for the G77 compiler).
       do i=1,mlabel
         collabel(i:i)=' '
       end do

      write(6,'(a)')
     $' vsf [t], vertical column [v], los column [l],',
     $' continuum [c], tilt [m], frequency shift [f],', 
     &'  solar-gas frequency shift [s], rms [r]: '
      read(5,*) ans

c  Find the number of windows/columns (NCOL)
      open(lunm,file='multiggg.bat',status='old')
      do icol=1,mcol  
         read(lunm,'(a)',end=99) tabel
      end do  ! icol=1,mcol     !  main loop (over windows)
      read(lunm,*,end=99) tabel
      write(6,*) 'Increase parameter mcol'
 99   close(lunm)
      ncol=icol-1

c  Read in the retrieved absorber amounts (YOBS+-YERR)
      open(lunm,file='multiggg.bat',status='old')
      open(lunv,file='collate.rpt',status='unknown')
      write(lunv,*)' iwin     window                    Nrow     Npp
     &NIT  %Conv  RMS_Min  RMS_Mean  RMS_Max'
      do icol=1,ncol     !  main loop (over windows)
        npp=0
        trms=0.0d0
        totnit=0
        ntc=0
        jval=icol-ncol
135     read(lunm,'(a)') tabel
        if(tabel(1:1).eq.':') go to 135
        colfile=tabel(index(tabel,'<')+1:index(tabel,'.ggg'))//'col'
        write(*,*) colfile
        open(lunc,file=colfile,status='old') ! .col file
        idot=index(colfile,'.')
        i1=index(colfile(:idot),'_')
        if(i1.eq.0) i1=index(colfile(:idot),'^')   !  new format
        if(i1.eq.0) i1=index(colfile(:idot),'-')   !  old format
c        read(colfile(i1+1:idot-1),*) fcen
        collabel=collabel(:lnbc(collabel)+2)//colfile(:idot-1)//' '//
     $  colfile(:idot-1)//'_error'
        window(icol)=colfile(:idot-1)
c
c  Read header lines of .col file and locate column containing "OVC_gas".
c  in order to read data from appropriate target gas.
        read(lunc,*) nlhead
        read(lunc,'(a)') gfit_version
        read(lunc,'(a)') gsetup_version
        do k=4,nlhead-1
        read(lunc,'(a)')colabel
        if(index(colabel,'runlogs').gt.0) runlog=colabel
        end do
        read(colabel,*) fcen, width, mit
        read(lunc,'(a)')colabel
        ktg=1+index(colabel,' OVC_'//colfile(:i1-1))
        if ( ktg .gt. 1) then
          call substr(colabel(:ktg-1),cdum,1,nss)
          ktg=(nss-4)/4
        endif
        iyrwas=-99999
        doywas=-99999
        zpdwas=-99999.9d0
        irow=0
        rmin=9999999.0
        rmax=0.0
        spname_rl=' '

        lr=lnbc(runlog)
        if(runlog(lr-2:lr-2).eq.'o') then
           max_delta_t=0.0004  ! 1.44s (ACE)
        else
           max_delta_t=0.0025  ! 9.0s (ground-based TCCON)
        endif
c
c  Read auxilliary measurements from runlog
        open(lunr,file=runlog, status='old')   !DG000906
        read(lunr,*)
        r8was=-9999999.9d0
        do while (irow.le.mrow)   !  Loop over runlog records with different times
           specname_col='='
           read(lunc,'(a)',end=24) col_string
c           write(*,'(a)')col_string(:lnbc(col_string))
           l2=fbc(col_string(2:))+1       ! First space following spectrum name
           l3=fnbc(col_string(l2:))+l2-1  ! First character of NIT
           l4=fbc(col_string(l3:))+l3-1   ! First space following NIT
c           write(*,*)l2,l3,l4

           write(csformat,'(a,i2.2,a)')'(1x,a',l4-5,
     &     ',i3,f6.3,3f5.1,f6.3,f7.4,f8.3,15(f7.3,e11.4,f9.4,e8.1))'
c           write(csformat,'(a,i2.2,a)')'(1x,a',l4-5,
c     &     ',i3,f6.3,3f5.1,f6.3,f7.4,f7.3,9(f7.3,e10.4,f9.4,e8.1))'

c           write(*,*) csformat
           read(col_string,csformat) specname_col,nit,cl,tilt,fqshift,
     &     sg,zlo,rmsfit,zmin,(airmass,ovcol,vsf,vsf_err,jtg=1,ktg)
           totnit=totnit+nit
           if(nit.lt.mit) ntc=ntc+1  ! Number of Times Converged
           if(rmsfit.le.0.0) then
              write(*,*) colfile,irow
              write(*,*)specname_col,nit,cl,tilt,fqshift,sg,zlo,rmsfit
              stop 'rmsfit <= 0'
           endif
           if(rmsfit.gt.rmax) rmax=rmsfit
           if(rmsfit.lt.rmin) rmin=rmsfit
24         continue

           do while(specname_col.ne.spname_rl)
              call read_runlog(lunr,col1,spname_rl,iyr,doy,zpdtim,
     &        obslat,obslon,zobs,asza,zenoff,opd,
     &        fovi,fovo,amal,ifirst,ilast,
     &        graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     &        tout,pout,hout,lasf,wavtkr,sia,sis,aipl,istat)
c              if(istat.ne.0) write(*,*) 'istat=',istat
c              write(*,*) '...called read_runlog.'
              if(istat.ne.0) go to 14     ! Exit Loop  irow=1,mrow
              lrl=lnbc(spname_rl)
              delta_t=zpdtim-zpdwas
c              write(*,*)spname_rl(:lrl),delta_t, max_delta_t
              if(iyr.ne.iyrwas .or. doy.ne.doywas .or.
     &        dabs(delta_t).ge.max_delta_t ) then  ! New observation time
                 if(spname_rl(lrl-2:lrl).ne.spname_rlwas(lrl-2:lrl))then
                    irow=irow+1
                    jval=jval+ncol   !  jval=icol+ncol*(irow-1)
                    yobs(jval)=ymiss
                    yerr(jval)=ymiss
                 endif
c              else
c                 write(*,*)'Warning: Delta_t < max_delta_t'
c                 write(*,*)spname_rl(:lrl),delta_t, max_delta_t
              endif
              iyrwas=iyr
              doywas=doy
              zpdwas=zpdtim
              spname_rlwas=spname_rl
           end do     ! while(specname_col.ne.spname_rl)
        
           sign(irow)=col1
           spectrum(irow)=spname_rl
           r8year=iyr+(doy+zpdtim/24.0d0)/365.25d0
           r8ydiff=r8year-r8was
           if( r8ydiff .lt. -0.00000001d0) write(*,*)
     &   'Warning:  Negative time step (runlog unsorted?) ',
     *    spname_rl,iyr,doy,zpdtim,r8was,r8year
           r8was=r8year
           year(irow)=r8year
c           yaux(1,irow)=r8year
           yaux(2,irow)=doy+zpdtim/24.0d0
           yaux(3,irow)=zpdtim
           yaux(4,irow)=irow
           yaux(5,irow)=obslat
           yaux(6,irow)=obslon
           yaux(7,irow)=zobs
           yaux(8,irow)=pout
           yaux(9,irow)=tout
           yaux(10,irow)=hout
           yaux(11,irow)=pins
           yaux(12,irow)=tins
           yaux(13,irow)=asza+zenoff
           yaux(14,irow)=graw
           yaux(15,irow)=zmin
           yaux(16,irow)=opd
           yaux(17,irow)=sqrt(fovi**2+amal**2)
           yaux(18,irow)=sia
           yaux(19,irow)=sis
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
           elseif(ans.eq.'f') then
               yobs(jval)=1000*fqshift/fcen
               yerr(jval)=vsf_err
c             yerr(jval)=vsf_err*abs(yobs(icol,irow))
           elseif(ans.eq.'s') then
               yobs(jval)=sg
               yerr(jval)=vsf_err
           elseif(ans.eq.'m') then
               yobs(jval)=tilt
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
        end do   !  irow.lt.mrow  !  Loop over spectra
        stop 'Increase parameter MROW'
14      nrow=irow
        if(ncol*nrow.gt.mval) then
           write(*,*)ncol,nrow
           write(*,*)'Increase parameter MVAL to ',nrow*ncol
           stop 'collate'
        endif
        close(lunr)
        close(lunc)
        write(lunv,'(i3,a30,2i8,2f7.2,3f9.4)')icol,'  '//colfile,
     &  npp, nrow,
     &  float(totnit)/npp,100*float(ntc)/npp,rmin,1/sqrt(trms/npp),rmax
      end do  ! icol=1,mcol     !  main loop (over windows)
      close(lunm)
      close(lunv)
c====================================================================
      do k=lr,1,-1
         if(ichar(runlog(k:k)) .eq. 92) go to 101  ! backslash
         if(ichar(runlog(k:k)) .eq. 47) go to 101  ! forward slash
      end do
101   outfile=runlog(k+1:lr-3)//ans//'sw'
c==================================================================
c  Write out all analyzed abundances to the .?sw disk file
      open(lunw,file=outfile,status='unknown')
      write(lunw,'(i2,i4,i6)') 6,nauxcol+2*ncol,nrow
      write(lunw,'(a)') collate_version
      write(lunw,'(a)') gfit_version
      write(lunw,'(a)') gsetup_version
c      write(lunw,'(a8,<nauxcol+2*ncol>(1pe12.4))') 'MISSING:',
      write(lunw,'(a8,1016(1pe12.4))') 'MISSING:',
     $(ymiss,j=1,nauxcol+2*ncol)
      write(lunw,'(a)')'  year  day  hour  run  lat  long  zobs  pout'//
     $' tout  hout  pins  tins  asza  graw  zmin  opd  fovi  sia  sis'//
     & collabel(:lnbc(collabel)+1)
c  Note that the SIGN array and the following IF statement are merely
c  to support both the new and the old runlog formats.
      nfound=0
      nmiss=0
      jval=0
      do irow=1,nrow
          do k=1,ncol
             jval=jval+1     !    jval=k+ncol*(irow-1)
             if( yobs(jval).eq.ymiss .and. yerr(jval).eq.ymiss ) then
             write(*,*)'missing', window(k),'  '//spectrum(irow)
                nmiss=nmiss+1
             else
                nfound=nfound+1
             endif
          end do
          write(lunw,75) sign(irow),year(irow),
     &    (yaux(k,irow),k=2,nauxcol),
     $    (yobs(k+ncol*(irow-1)),yerr(k+ncol*(irow-1)),k=1,ncol)
75        format(a1,f13.8,18f13.5,800(1pe12.4))
      end do  !  irow=1,nrow
      close(lunw)
c====================================================================
      write(6,*)
      write(6,*) outfile//' contains:'
      write(6,'(i7,a)') nauxcol,' auxiliary columns'
      write(6,'(i7,a)') ncol,' data columns (each value + error)'
      write(6,'(i7,a)') nrow,' data rows'
      write(6,'(i7,a,f5.1,a1)') nfound,' found values   = ',
     & 100*float(nfound)/(nfound+nmiss),'%'
      write(6,'(i7,a,f5.1,a1)') nmiss,' missing values = ',
     & 100*float(nmiss)/(nfound+nmiss),'%'
c====================================================================
      stop
      end
