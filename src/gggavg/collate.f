c  Program COLLATE
c
c  GGGAVG reads the .col output files (gas_1234.runlog.col) produced by GFIT,
c  collates them with auxiliary measurements from the runlog, and then writes
c  the information to a spreadsheet-format XYPLOT-readable output file named
c  runlog.tsw, runlog.csw, etc.
c
c  INPUT FILES:
c     runlog.xrl
c     gas_1234.runlog.col
c     multiggg.bat     batch file containing names of .col files
c
c OUTPUT FILES:
c     runlog.xsw   Spreadsheet of individual window values for each spectrum
c     runlog.xav   Spreadsheet of averaged (over each gas) values for each spectrum
c     wincomp.xav  Statistical information, obtained during averaging, about
c                  the consistency of different windows of the same gas.

      implicit none
      integer i1,idot,irow,j,jtg,k,ktg,lnbc,fnbc,fbc,
     & npp,lrl,l2,l3,l4,totnit,ntc,
     &lr,lunm,lunr,lunc,luno,lunt,lunv,lunw,mcol,ncol,icol,
     &mrow,nauxcol,nlhead,nmiss_av,nmiss_sw,iyrwas,doywas,
     &nrow,nss,num,nit,mit,i,iyr,doy,istat,nfound_sw,nfound_av
      parameter (lunm=11)
      parameter (lunr=12)
      parameter (lunc=13)
      parameter (lunt=14)
      parameter (lunw=15)
      parameter (luno=16)
      parameter (lunv=17)
      parameter (mcol=600)     ! Total number of columns/windows
      parameter (mrow=160000)  ! Max number of output records/spectra
      parameter (nauxcol=19)   ! Number of auxiliary parameters/columns

      character ans*1,apf*2,cdum*20,colabel*500,
     & gfit_version*80,gsetup_version*80,col_string*500,
     & csformat*80,cl_sw*16000,outfile_sw*20,outfile_av*20,col1*1,
     & spname_rl*21,specname_col*21,runlog*72,sign(mrow)*1,
     & spectrum(mrow)*20,tabel*80,
     & colfile*40,collate_version*44,window(mcol)*10,spname_rlwas*21

      real*8 airmass,asza,cl,tilt,zlo,error_term,
     & fcen,width,zobs,rmin,rmax,
     & fqshift,graw,obslat,obslon,opd,ovcol,rmsfit,
     & tew,totcon,
     & toterr,ymiss,yaux(nauxcol,mrow),yerr(mcol,mrow),
     & r8was,r8year,r8ydiff,
     & yobs(mcol,mrow),weight,twt,tsc,trms,lasf,wavtkr,sia,sis,aipl
      parameter (ymiss=9.9999d+29)

      real*8 zpdtim,tout,pout,hout,tins,pins,hins,fovi,fovo,amal,
     & snr,zenoff,zoff,sg,zpdwas,max_delta_t,delta_t,zmin
      integer bytepw,ifirst,ilast,possp

      collate_version=' COLLATE    Version 1.0.0   14-Sep-2008   GCT'
      write(6,*) collate_version

c  Initialize character arrays (Necessary for the G77 compiler).
       gas(1)='........'
       do i=1,16000
         cl_sw(i:i)=' '
       end do

      write(6,'(a)')
     $' totcon [t], vertical column [v], los column [l],',
     $' continuum [c], tilt [m], frequency shift [f],', 
     &'  solar-gas frequency shift [s], rms [r]: '
      read(5,*) ans

c  Read in the retrieved absorber amounts (YOBS+-YERR)
      cl_sw=' '
      open(lunm,file='multiggg.bat',status='old')
      open(lunv,file='collate.rpt',status='unknown')
      write(lunv,*)' iwin     window                    Nrow     Npp
     &NIT  %Conv  RMS_Min  RMS_Mean  RMS_Max'
      do icol=1,mcol     !  main loop (over windows)
        npp=0
        trms=0.0d0
        totnit=0
        ntc=0
135     read(lunm,'(a)',end=99) tabel
        if(tabel(1:1).eq.':') go to 135
        colfile=tabel(index(tabel,'<')+1:index(tabel,'.ggg'))//'col'
        write(*,*) colfile
        open(lunc,file=colfile,status='old') ! .col file
        idot=index(colfile,'.')
        i1=index(colfile(:idot),'_')
        if(i1.eq.0) i1=index(colfile(:idot),'^')   !  new format
        if(i1.eq.0) i1=index(colfile(:idot),'-')   !  old format
c        read(colfile(i1+1:idot-1),*) fcen
        cl_sw=cl_sw(:lnbc(cl_sw)+2)//colfile(:idot-1)//' '//
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
     &     sg,zlo,rmsfit,zmin,(airmass,ovcol,totcon,toterr,jtg=1,ktg)
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
                    yobs(icol,irow)=ymiss
                    yerr(icol,irow)=ymiss
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
     &   'Warning:  Negative time step (runlog sorted?)',
     *    spname_rl,iyr,doy,zpdtim,r8was,r8year
           r8was=r8year
           yaux(1,irow)=r8year
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
               yobs(icol,irow)=totcon
               yerr(icol,irow)=toterr
           elseif(ans.eq.'v') then
               yobs(icol,irow)=totcon*ovcol
               yerr(icol,irow)=toterr*ovcol
           elseif(ans.eq.'l') then
               yobs(icol,irow)=totcon*ovcol*airmass
               yerr(icol,irow)=toterr*ovcol*airmass
           elseif(ans.eq.'f') then
               yobs(icol,irow)=1000*fqshift/fcen
               yerr(icol,irow)=toterr
c             yerr(icol,irow)=toterr*abs(yobs(icol,irow))
           elseif(ans.eq.'s') then
               yobs(icol,irow)=sg
               yerr(icol,irow)=toterr
           elseif(ans.eq.'m') then
               yobs(icol,irow)=tilt
               yerr(icol,irow)=toterr
           elseif(ans.eq.'r') then
               yobs(icol,irow)=rmsfit
               yerr(icol,irow)=toterr*rmsfit
           elseif(ans.eq.'c') then
               yobs(icol,irow)=cl
               yerr(icol,irow)=cl*rmsfit
           else
               stop 'unknown option'
           endif
           npp=npp+1
           trms=trms+1.0d0/rmsfit**2
        end do   !  irow.lt.mrow  !  Loop over spectra
        stop 'Increase parameter MROW'
14      nrow=irow
        close(lunr)
        close(lunc)
        write(lunv,'(i3,a30,2i8,2f7.2,3f9.4)')icol,'  '//colfile,
     &  npp, nrow,
     &  float(totnit)/npp,100*float(ntc)/npp,rmin,1/sqrt(trms/npp),rmax
      end do  ! icol=1,mcol     !  main loop (over windows)
      read(lunm,*,end=99) tabel
      write(6,*) 'Increase parameter mcol'
 99   close(lunm)
      close(lunv)
      ncol=icol-1
c====================================================================
      do k=lr,1,-1
         if(ichar(runlog(k:k)) .eq. 92) go to 101  ! backslash
         if(ichar(runlog(k:k)) .eq. 47) go to 101  ! forward slash
      end do
101   outfile_sw=runlog(k+1:lr-3)//ans//'sw'
c==================================================================
c  Write out all analyzed abundances to the .?sw disk file
      open(lunw,file=outfile_sw,status='unknown')
      write(lunw,'(2i3)') 6,nauxcol+2*ncol
      write(lunw,'(a)') collate_version
      write(lunw,'(a)') gfit_version
      write(lunw,'(a)') gsetup_version
c      write(lunw,'(a8,<nauxcol+2*ncol>(1pe12.4))') 'MISSING:',
      write(lunw,'(a8,1016(1pe12.4))') 'MISSING:',
     $(ymiss,j=1,nauxcol+2*ncol)
      write(lunw,'(a)')'  year  day  hour  run  lat  long  zobs  pout'//
     $' tout  hout  pins  tins  asza  graw  zmin  opd  fovi  sia  sis'//
     & cl_sw(:lnbc(cl_sw)+1)
c  Note that the SIGN array and the following IF statement are merely
c  to support both the new and the old runlog formats.
      nfound_sw=0
      nmiss_sw=0
      do irow=1,nrow
          do k=1,ncol
             if( yobs(k,irow).eq.ymiss .and.
     $       yerr(k,irow).eq.ymiss ) then
             write(*,*)'missing',window(k),'  '//spectrum(irow-1)
                nmiss_sw=nmiss_sw+1
             else
                nfound_sw=nfound_sw+1
             endif
          end do
          write(lunw,75) sign(irow),
     &    (yaux(k,irow),k=1,nauxcol),
     $    (yobs(k,irow),yerr(k,irow),k=1,ncol)
      end do  !  irow=1,nrow
      close(lunw)
c====================================================================
      write(6,*)
      write(6,*) outfile_sw//' contains:'
      write(6,'(i7,a)') nauxcol,' auxiliary columns'
      write(6,'(i7,a)') ncol,' data columns (each value + error)'
      write(6,'(i7,a)') nrow,' data rows'
      write(6,'(i7,a,f5.1,a1)') nfound_sw,' found values   = ',
     & 100*float(nfound_sw)/(nfound_sw+nmiss_sw),'%'
      write(6,'(i7,a,f5.1,a1)') nmiss_sw,' missing values = ',
     & 100*float(nmiss_sw)/(nfound_sw+nmiss_sw),'%'
c====================================================================
      stop
      end
