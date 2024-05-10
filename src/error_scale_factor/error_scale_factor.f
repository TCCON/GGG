c  Program: error_scale_factor.f
c
c  Purpose: To determine the factor by which the error bars
c  returned by GFIT/COLLATE/AVERAGE exceed the true precision
c  of the measurements in the selected runlog.vav.ada file
c
c  Input Files:
c       runlog.vav.ada 
c
c  Output:
c       Written to screen
c       
c  Due to errors/omissions in the spectroscopy and perhaps also
c  in the modeling of the instrument response (e.g. continuum, ILS),
c  the spectral fitting residuals produced by GFIT tend to be dominated
c  by systematic artifacts, not random noise. These artifacts are generally
c  the same from spectrum to spectrum, so an error calculated from the
c  rms spectral fit will generally overestimate the true precision.
c
c  This program computes the difference in the column amounts between
c  sucessive spectra and compares this with the quadrature sum of their 
c  errors. The ratio should average out to 1.0 if the error bars are
c  consistent with the scatter of the measured values. The average value
c  of this ratio is calculated for each gas, the weights being given by
c  a Lorentzian with a 5 minute half-width.
c
c  Uncertainties are also computed for the daily error scale factors.
c  These are based on the sum of the weights for each day. For a day
c  with a single pair of observation very close in time, the ESF
c  uncertainty is assumed to be 100%. The uncertainty reduces as the
c  number of comparable observations increases, The uncertainty
c  increases as the comparable observations are less simultaneous.
c  The fractional uncertainty is the same for every gas, because it
c  depends only on the number and the timing of the observations.
c
      implicit none
      include "../gfit/ggg_int_params.f"
      include "../comn/postproc_params.f"

      integer*4 lunr,ncoml,ncol,mcol,kcol,j,icount,
     & lunw_d,lunw_g,
     & lnbc,naux,nrow,li,k,kspflag,kday,kwas,istat,idum
      parameter (lunr=14,lunw_d=15,lunw_g=16,mcol=150)
      character header*2400,headarr(mcol)*12,specname*(nchar),
     & inputfile*64,alloutputfile*80,dailyoutputfile*80,version*62,
     & cdum*(mcharhead), headdum(maddln)*(mcharhead)
      real*8 yrow(mcol),ywas(mcol),dt,wt,
     & tdwt,totd(mcol),yy,dtiny,
     & tgwt,totg(mcol),rdum

      idum=mfilepath ! Avoid compiler warning (unused parameter)
      idum=mauxcol  ! Avoid compiler warning (unused parameter)
      idum=mcolvav  ! Avoid compiler warning (unused parameter)
      idum=mgas     ! Avoid compiler warning (unused parameter)
      idum=mlev     ! Avoid compiler warning (unused parameter)
      idum=mrow_qc  ! Avoid compiler warning (unused parameter)
      idum=mspeci   ! Avoid compiler warning (unused parameter)
      idum=mvmode   ! Avoid compiler warning (unused parameter)
      idum=ncell    ! Avoid compiler warning (unused parameter)
      idum=nchar    ! Avoid compiler warning (unused parameter)

      version=
     &' error_scale_factor        Version 1.23a  2020-03-12   GCT,JLL'
      write(*,*)version
      dtiny=0.000001d0

      kspflag=0
      icount=1  ! avoid compiler warnings (may be uninitialized)
      kwas=0    ! avoid compiler warnings (may be uninitialized)
      tdwt=0.0  ! avoid compiler warnings (may be uninitialized)
      tgwt=0.0

      if (iargc() == 0) then
         write(*,*)'Enter name of input file (e.g. paIn_1.0lm.vav.ada):'
         read(*,'(a)') inputfile
      elseif (iargc() == 1) then
         call getarg(1, inputfile)
      else
         stop 'Usage: $gggpath/bin/error_scale_factor adafile'
      endif

      li=lnbc(inputfile)
      open(lunr,file=inputfile, status='old')

c  Read the header of the .ada file. Most of the information
c  is not needed so we store it in dummy variables.
      call read_postproc_header(lunr, ncoml, ncol, nrow, naux,
     & rdum, cdum, headdum, idum)
      read(lunr,'(a)') header
      call substr(header,headarr,mcol,kcol)
      call lowercase(header)
      
      if(ncol.gt.mcol) stop 'increase mcol'
      do k=1,ncol
         totd(k)=0.0
         totg(k)=0.0
      end do
      if (index(header,'spectrum') .gt. 0) kspflag=1
      if(kcol.ne.ncol ) stop 'ncol/kcol mismatch'

      alloutputfile=inputfile(:li)//'.all_error.out'
      open(lunw_g,file=alloutputfile,status='unknown')
      write(lunw_g,*) 2,1+(ncol-naux)/2
      write(lunw_g,'(6x,64a12)')' Year       ',
     & (headarr(k),k=naux+1,ncol-1,2)

      dailyoutputfile=inputfile(:li)//'.daily_error.out'
      open(lunw_d,file=dailyoutputfile,status='unknown')
      write(lunw_d,*) 2,3+(ncol-naux)
      write(lunw_d,'(6x,2a12,a5,40(a12,a20))')' Year       ',
     & ' DOY       ',' N   ',
     & (headarr(k),headarr(k)(:lnbc(headarr(k)))//'_error  ',
     & k=naux+1,ncol-1,2)

c  Read each day of data into memory and divide XGas values by the
c  appropriate correction factors.
      dt=0.00001  ! years (=5 minutes)
      if (kspflag .eq. 1) then
         read(lunr,*) specname, (ywas(j),j=1+kspflag,ncol)
      else
         read(lunr,*) (ywas(j),j=1,ncol)
      endif
      kwas = nint(ywas(2+kspflag)-ywas(3+kspflag)/24)

      istat=0
      do while (istat.eq.0)
         if (kspflag .eq. 1) then
            read(lunr,*,iostat=istat)specname,(yrow(j),j=1+kspflag,ncol)
         else
            read(lunr,*,iostat=istat) (yrow(j),j=1,ncol)
         endif
         kday = nint(yrow(2+kspflag)-yrow(3+kspflag)/24)
c         write(*,*) specname(:22),
c     &   yrow(2+kspflag),yrow(3+kspflag),kday,kwas

         if(kday+istat.eq.kwas) then
            icount=icount+1
c  Lorentzian weights based on time difference (first column).
            wt=1/(1+((yrow(1+kspflag)-ywas(1+kspflag))/dt)**2)  ! Lorentzian Weights
            tdwt=tdwt+wt
            do k=naux+1,ncol-1,2
               yy=(yrow(k)-ywas(k))**2/(yrow(k+1)**2+ywas(k+1)**2)
               totd(k)=totd(k)+wt*yy
            end do
         else
c  Output daily results. If there's only one observation, 
c  skip because the ESF will be indeterminate (NaN)
c  Don't let ESF fall below SQRT(0.0001) = 0.01
c  for days when the measured columns are all identical
c            if(icount.ge.2) then
            write(lunw_d,'(2f12.6,1x,i4,70f12.6)')
     &      ywas(1+kspflag),
     &      ywas(2+kspflag),icount,
     &      (sqrt(0.0001+totd(k)/tdwt),
     &      sqrt(0.01+totd(k)/tdwt)/sqrt(tdwt),
     &      k=naux+1,ncol-1,2)
c            endif

c  Augment global totals
            tgwt=tgwt+tdwt
            do k=1,ncol
               totg(k)=totg(k)+totd(k)
            end do

c  Reset daily totals
            icount=1
            tdwt=dtiny
            do k=1,ncol
               totd(k)=dtiny
            end do

         endif     !  if(kday.eq.kwas) then

c  Copy current data to YWAS
         kwas=kday
         do k=1+kspflag,ncol
            ywas(k)=yrow(k)
         end do

      end do         ! do while (istat.eq.0)
      close(lunr)
      close(lunw_d)

c  Output global results.
      write(lunw_g,'(64f12.6)')ywas(1+kspflag),
     & (min(9999.0,sqrt(totg(k)/tgwt)),k=naux+1,ncol-1,2)
      close(lunw_g)

c  Output results to the screen.
      do k=naux+1,ncol-1,2
         if(headarr(k).ne.'luft') then ! avoid writing out the error
c     scale factor for column luft which isn't terribly meaningful
            write(*,'(i4,2x,a,f9.3)') k,headarr(k),sqrt(totg(k)/tgwt)
c            write(55,'(i4,2x,a,f9.3)') k,headarr(k),sqrt(totg(k)/tgwt)
         endif
      end do

      stop
      end
