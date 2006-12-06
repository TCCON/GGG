c  GGGAVG reads the output files (gas_1234.runlog.col) produced by GFIT,
c  collates them with auxiliary measurements from the runlog, and then
c  outputs the collated information to a spreadsheet-format file which is
c  readable by the IDL program AMESPLOT. It further uses subroutine TFIT
c  to average sucessive windows having the same target GAS, and also writes
c  these averaged values to an AMESPLOT readable spreadsheet file.
c
c  INPUT FILES:
c     runlog.xxx
c     gas_1234.runlog.col
c     multiggg     batch file containing names of .col files
c
c OUTPUT FILES:
c     runlog.xsw   Spreadsheet of individual window values for each spectrum
c     runlog.xav   Spreadsheet of averaged (over each gas) values for each spectrum
c     wincomp.xav  Statistical information, obtained during averaging, about
c                  the consistency of different windows of the same gas.

      implicit none
      integer i1,iavg,idot,irow,j,jtg,k,kavg,ktg,lnbc,
     $lr,lunm,lunr,lunc,luno,lunt,lunw,mavg,mcol,ncol,icol,
     $mrow,nauxcol,navg,nlhead,nn,nmissav,nmissing,nmisssw,
     $nrow,nss,num,nit,i,iyr,doy,istat
      parameter (lunm=11)
      parameter (lunr=12)
      parameter (lunc=13)
      parameter (lunt=14)
      parameter (lunw=15)
      parameter (luno=16)
      parameter (mavg=50)
      parameter (mcol=600)   ! number of windows
      parameter (mrow=20000)  ! number of spectra
      parameter (nauxcol=15)

      character ans*1,apf*2,cl_av*1200,cdum*20,colabel*500,
     $gas(0:mavg)*8,gfitversion*80,gsetupversion*80,
     $cl_sw*16000,outfile_sw*20, outfile_av*20,
     $runlab*20,runlog*72,sign(mrow),spectrum(mrow)*20,tabel*80,
     $colfile*40,version*44,window(mcol)*10

      real*8 airmass,asza,cl,cle,tilt,zlo,error_term,fcen,zobs,
     $fqshift,graw,obslat,obslon,opd,ovcol,rmsfit,scale(mcol),
     $sew(mcol),sigma_scale(mcol),sigma_vmr(mrow),tew,totcon,
     $toterr,valmiss,vmr(mrow),yaux(nauxcol,mrow),yerr(mcol,mrow),
     $yobs(mcol,mrow),weight,twt,tsc,trms,lasf,wavtkr
      parameter (valmiss=9.9999d+29)

      real*8 zpdtim,tout,pout,hout,tins,pins,hins,fovi,fovo,amal,
     & snr,zenoff,zoff,sg
      integer avindx(mavg+1),bytepw,ifirst,ilast,possp

      version=' GGGAVG    Version 2.4.0   11-Mar-2004   GCT'
      write(6,*) version
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c     Platform specification: DG000909
c      call getenv('LOGNAME',user)
c      if(user.ne.'        ')then
c         platform=0               !0=Sun, 1=PC-Linux, 2=PC-Win32
c         dl='/'
c         root='/ggg/'
c         batchfile='multiggg'
c      else
c         platform=2               !0=Sun, 1=PC-Linux, 2=PC-Win32
c         dl=char(92)  ! backslash ('\')
c         root='c:'//dl//'ggg'//dl
c         user='PC-Win'
c         batchfile='multiggg.bat'
c      endif
c      lrt=lnbc(root)       !Length of root
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c  Initialize character arrays (Necessary for the G77 compiler).
       gas(1)='........'
       do i=1,1200
         cl_av(i:i)=' '
       end do
       do i=1,9600
         cl_sw(i:i)=' '
       end do

      write(6,'(a,$)')
     $' totcon [t], vertical column [v], los column [l],',
     $' continuum [c], tilt [m], frequency shift [f], rms [r]: '
      read(5,*) ans

c  Read in the retrieved absorber amounts (YOBS+-YERR)
      navg=0
      nmisssw=0
      cl_av=' '
      cl_sw=' '
      open(lunm,file='multiggg.bat',status='old')	!DG000906
      do icol=1,mcol     !  main loop (over windows)
        trms=0.0d0
 80     read(lunm,'(a)',end=99) tabel
        colfile=tabel(index(tabel,'<')+1:index(tabel,'.ggg'))//'col'
        open(lunc,file=colfile,status='old',err=80) ! .col file
        idot=index(colfile,'.')
        i1=index(colfile(:idot),'_')
        if(i1.eq.0) i1=index(colfile(:idot),'^')   !  new format
        if(i1.eq.0) i1=index(colfile(:idot),'-')   !  old format
        read(colfile(i1+1:idot-1),*) fcen
        if(colfile(:i1-1).ne.gas(navg)) then ! ITS A NEW GAS
          if(navg.lt.mavg) then
            navg=navg+1
            gas(navg)=colfile(:i1-1)
            avindx(navg)=icol
            write(6,*) navg,avindx(navg)
            cl_av=cl_av(:lnbc(cl_av)+2)//colfile(:i1-1)//' '//
     $      colfile(:i1-1)//'_error'
          else
            write(6,*) 'Increase parameter MAVG >',navg
            stop
          endif
        endif
        cl_sw=cl_sw(:lnbc(cl_sw)+2)//colfile(:idot-1)//' '//
     $  colfile(:idot-1)//'_error'
        window(icol)=colfile(:idot-1)
c
c  Read header lines of .col file and locate column containing "OVC_gas".
c  in order to read data from appropriate target gas.
        read(lunc,*) nlhead,nn
        read(lunc,'(a)') gfitversion
        read(lunc,'(a)') gsetupversion
        do k=4,nlhead
        read(lunc,'(a)')colabel
        if(index(colabel,'runlogs').gt.0) runlog=colabel
        end do
        ktg=1+index(colabel,' OVC_'//colfile(:i1-1))
        if ( ktg .gt. 1) then
          call substr(colabel(:ktg-1),cdum,1,nss)
          write(*,*)'ktg=',ktg,nss
          ktg=(nss-4)/4
        endif
        nmissing=0
c
c  Read auxilliary measurements from runlog
        lr=lnbc(runlog)
        open(lunr,file=runlog, status='old')   !DG000906
        read(lunr,*)
        do irow=1,mrow  !  Loop over spectra
 12     call readrunlog(lunr,sign(irow),spectrum(irow),iyr,doy,
     &  zpdtim,obslat,obslon,zobs,asza,zenoff,opd,fovi,fovo,amal,
     &  ifirst,ilast,graw,possp,bytepw,zoff,snr,apf,tins,pins,hins,
     &  tout,pout,hout,lasf,wavtkr,istat)
        if(istat.ne.0) go to 14

c 12       read(lunr,'(a)',end=14) record
c          if(record(1:1).eq.':') go to 12
c          iflag=1
c          if(record(1:1).eq.' '.or.record(1:1).eq.'+'.or.
c     $    record(1:1).eq.'-') then
c            iflag=2
c            sign(irow)=record(1:1)
c          endif
c          if(runlog(lr-2:lr).eq.'oat') then  ! ATMOS runlogs are TAB delimited
c            read(record(iflag:),*) spectrum(irow),year,doy,zpdtim,
c     $      obslat,obslon,zobs,asza,zenoff,opd,fovi,fovo,amal,ifirst,
c     $      ilast,graw,possp,bytepw,zoff,zerr,snr,scalf,apf
c          else
c            read(record(iflag:),333) spectrum(irow),year,doy,zpdtim,
c     $      obslat,obslon,zobs,asza,zenoff,opd,fovi,fovo,amal,
c     $      ifirst,ilast,graw,possp,bytepw,zoff,snr,apf,tout,pout,hout
c          endif
c 333      format(a20,2x,2i4,f8.4,f7.3,3f8.3,f7.0,f7.2,3f6.0,
c     $    2i8,f14.11,i8,i3,2f5.0,1x,a2,f4.0,f8.2,f3.0)

c          if(zpdtim.lt.2.0) zpdtim=zpdtim+24.  ! Park Falls Fudge
          yaux(1,irow)=iyr+(doy+zpdtim/24.0d0)/365.25d0
          yaux(2,irow)=doy+zpdtim/24.0d0
          yaux(3,irow)=zpdtim
          yaux(4,irow)=irow
          yaux(5,irow)=obslat
          yaux(6,irow)=obslon
          yaux(7,irow)=zobs
          yaux(8,irow)=pout
          yaux(9,irow)=tout
          yaux(10,irow)=pins
          yaux(11,irow)=asza+zenoff
          yaux(13,irow)=graw
          yaux(14,irow)=opd
          yaux(15,irow)=sqrt(fovi**2+amal**2)
c
 44       read(lunc,93,end=95) runlab,nit,cl,tilt,fqshift,sg,zlo,rmsfit,
     $    yaux(12,irow),(airmass,ovcol,totcon,toterr,jtg=1,ktg)
 93       format(1x,a20,i2,f6.3,f5.1,f4.1,f5.1,f6.3,f7.4,f7.3,
     $    9(f7.3,e10.3,f9.4,e8.1))
c 93       format(1x,a20,i2,f6.3,2f5.1,f6.3,f7.4,f7.3,
c     $    9(f7.3,e10.3,f9.4,e8.1))
c 93       format(1x,a20,i3,f6.3,2f5.1,f6.3,f7.4,f7.3,
c     $    9(f7.3,e10.3,f9.4,e8.1))
          if(runlab(2:3).eq.'hg'.or.runlab(2:3).eq.'in')
     $    runlab(1:3)=spectrum(irow)(1:3)
 95       if(runlab.ne.spectrum(irow)) then
            yobs(icol,irow)=valmiss
            yerr(icol,irow)=valmiss
            backspace(lunc)
            nmissing=nmissing+1
          else
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
              yobs(icol,irow)=0.001d0*fqshift/fcen
              yerr(icol,irow)=toterr
c              yerr(icol,irow)=toterr*abs(yobs(icol,irow))
            elseif(ans.eq.'m') then
              yobs(icol,irow)=tilt
              yerr(icol,irow)=toterr
            elseif(ans.eq.'r') then
              yobs(icol,irow)=rmsfit
              yerr(icol,irow)=toterr*rmsfit
            elseif(ans.eq.'c') then
              yobs(icol,irow)=cl
              cle=0.0
              yerr(icol,irow)=cle
            else
              stop 'unknown option'
            endif
          endif    !  if(runlab.ne.spectrum(irow))
          trms=trms+1.0d0/rmsfit**2
        end do   !  irow=1,nrow  !  Loop over spectra
        stop 'Increase parameter MROW'
14      nrow=irow-1
        close(lunr)
        close(lunc)
        write(6,'(i4,a,2i6,f10.6)')
     &  icol,'  '//colfile, nrow-nmissing,nrow,1/sqrt(trms/nrow)
        nmisssw=nmisssw+nmissing
      end do  ! icol=1,mcol     !  main loop (over windows)
      read(lunm,*,end=99) tabel
      write(6,*) 'Increase parameter mcol'
 99   close(lunm)
      ncol=icol-1
      avindx(navg+1)=icol
c====================================================================
      do k=lr,1,-1
         if(ichar(runlog(k:k)) .eq. 92) go to 101  ! backslash
         if(ichar(runlog(k:k)) .eq. 47) go to 101  ! forward slash
      end do
101   outfile_sw=runlog(k+1:lr-3)//ans//'sw'
      outfile_av=runlog(k+1:lr-3)//ans//'av'
c  Write out all analyzed abundances to disk file
      write(6,*)
      write(6,*) outfile_sw//' contains:'
      write(6,'(i5,a)') nauxcol,' auxiliary columns'
      write(6,'(i5,a)') ncol,' data columns (each value + error)'
      write(6,'(i5,a)') nrow,' data rows'
      write(6,'(i5,a)') nmisssw,' missing values'
      open(lunw,file=outfile_sw,status='unknown')
      write(lunw,*) 6,nauxcol+2*ncol
      write(lunw,'(a)') version
      write(lunw,'(a)') gfitversion
      write(lunw,'(a)') gsetupversion
c      write(lunw,'(a8,<nauxcol+2*ncol>(1pe12.4))') 'MISSING:',
      write(lunw,'(a8,1014(1pe12.4))') 'MISSING:',
     $(valmiss,j=1,nauxcol+2*ncol)
c==================================================================
      write(lunw,'(a)')'  year  day  hour  run  lat  long  zobs  pobs'//
     $'  tobs  pins  asza  zmin  graw  opd  fovi'//cl_sw(:lnbc(cl_sw)+1)
c  Note that the SIGN array and the following IF statement are merely
c  to support both the new and the old runlog formats.
      do irow=1,nrow
          write(lunw,75) sign(irow),(yaux(k,irow),k=1,nauxcol),
     $    (yobs(k,irow),yerr(k,irow),k=1,ncol)
      end do  !  irow=1,nrow
      close(lunw)
c====================================================================
c  Compute average absorber amounts for gases having more than 1 window
      open(luno,file='outliers',status='unknown')
      open(lunt,file='wincomp.'//ans//'av',status='unknown')
      write(lunt,*) 2 4
      write(lunt,*)' Window   Mean_Col   Std_Err   Chi2/N'
      nmissav=0
      do iavg=1,navg
        kavg=avindx(iavg)
        num=avindx(iavg+1)-kavg
        if(num.eq.1) then
          do irow=1,nrow
            yobs(iavg,irow)=yobs(avindx(iavg),irow)
            yerr(iavg,irow)=yerr(avindx(iavg),irow)
            if( yobs(iavg,irow).eq.valmiss .and.
     $      yerr(iavg,irow).eq.valmiss ) nmissav=nmissav+1
          end do
        elseif(num.gt.1) then
          call t_fit(valmiss,nrow,mcol,num,yobs(kavg,1),yerr(kavg,1),
     $    vmr,sigma_vmr,scale(kavg),sigma_scale(kavg),sew(kavg),tew)
          write(6,*) gas(iavg),' Total Error Weight =',tew
          twt=0.0d0
          tsc=0.0d0
          do icol=kavg,kavg+num-1
            write(lunt,'(a,3f9.4)') window(icol),scale(icol),
     $      sigma_scale(icol),sew(icol)
            weight=1.0d0/sigma_scale(icol)**2
            twt=twt+weight
            tsc=tsc+weight*scale(icol)
          end do
           write(lunt,*)'=================================='
           write(lunt,'(a,3f9.4)') 'Means     ',tsc/twt,1/sqrt(twt),tew
           write(lunt,*) 

          do irow=1,nrow
            do icol=kavg,kavg+num-1
              if(yerr(icol,irow).ne.valmiss) then
                error_term=(yobs(icol,irow)-vmr(irow)*scale(icol))/
     $          (tew*yerr(icol,irow))
                if(abs(error_term).gt.5.) write(luno,611)
     $          error_term,spectrum(irow),window(icol)
 611            format(' Deviation =',f9.5,' sigma  for spectrum ',a20,
     $          ' in window  ',a10)
              endif
            end do
            yobs(iavg,irow)=vmr(irow)
            if(vmr(irow).eq.valmiss.and.sigma_vmr(irow).eq.valmiss) then
              yerr(iavg,irow)=valmiss
              nmissav=nmissav+1
            else
              yerr(iavg,irow)=dmax1(1.0d0,tew)*sigma_vmr(irow)
            endif
          end do  ! irow=1,nrow
        else     ! num < 1
          stop 'error in indexing'
        endif
      end do     !  iavg=1,navg
      close(lunt)
      close(luno)
c=========================================================================
c  Write out average abundances to disk file
      write(6,*)
      write(6,*) outfile_av//' contains:'
      write(6,'(i5,a)') nauxcol,' auxiliary columns'
      write(6,'(i5,a)') navg,' data columns (each value + error)'
      write(6,'(i5,a)') nrow,' data rows'
      write(6,'(i5,a)') nmissav,' missing values'
      open(lunw,file=outfile_av,status='unknown')
      write(lunw,*) 6,nauxcol+2*navg
      write(lunw,'(a)') version
      write(lunw,'(a)') gfitversion
      write(lunw,'(a)') gsetupversion
c      write(lunw,'(a8,<nauxcol+2*navg>(1pe12.4))') 'MISSING:',
      write(lunw,'(a8,1014(1pe12.4))') 'MISSING:',
     $(valmiss,j=1,nauxcol+2*navg)
      write(lunw,'(a)')'  year  day  hour  run  lat  long  zobs  pobs'//
     $'  tobs  pins  asza  zmin  graw  opd  fovi'//cl_av(:lnbc(cl_av)+1)
      do  irow=1,nrow
          write(lunw,75) sign(irow),(yaux(k,irow),k=1,nauxcol),
     $    (yobs(iavg,irow),yerr(iavg,irow),iavg=1,navg)
 75       format(1a,f11.6,13f12.5,800(1pe12.4))
      end do
      close(lunw)
      stop
      end
