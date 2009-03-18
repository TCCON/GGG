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
c  Due to errors/ommissions in the spectroscopy and perhaps also
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
      implicit none
      integer*4 lunr,luns,ncoml,ncol,mcol,kcol,j,
     & lnbc,irow,naux,nrow,li,k, nchar
      parameter (lunr=14,luns=15,mcol=150)
      character header*800,headarr(mcol)*20,gggdir*80,specname*35,
     & inputfile*40,version*62
      real*8 yrow(mcol),ywas(mcol),dt,wt,twt,tot(mcol)

      version=
     &' error_scale_factor           Version 1.1.2   2009-03-03   GCT'

      nchar=0
      call getenv('GGGPATH',gggdir)

      write(*,*)'Enter name of input file (e.g. paIn_1.0lm.vav):'
      read(*,'(a)') inputfile
      li=lnbc(inputfile)
      open(lunr,file=inputfile, status='old')

c  Read the header of the .ada file and figure out the
      read(lunr,'(i2,i4,i7,i4)') ncoml,ncol,nrow,naux
      if(ncol.gt.mcol) stop 'increase mcol'
      do j=2,ncoml
         read(lunr,'(a)') header
      end do
      call substr(header,headarr,mcol,kcol)
      if (index(header,'Spectrum') .gt. 0) nchar=1
      if(kcol.ne.ncol ) stop 'ncol/kcol mismatch'

c  Read each day of data into memory and divide XGas values by the
c  appropriate correction factors.
      dt=0.00001  ! years (=5 minutes)
      twt=0.0
      do k=1,ncol
         tot(k)=0.0
      end do
      if (nchar .eq. 1) then
         read(lunr,*) specname, (ywas(j),j=1+nchar,ncol)
      else
         read(lunr,*) (ywas(j),j=1,ncol)
      endif
      do irow=1,9999999
         if (nchar .eq. 1) then
            read(lunr,*,end=99) specname, (yrow(j),j=1+nchar,ncol)
         else
            read(lunr,*,end=99) (yrow(j),j=1,ncol)
         endif

         wt=1/(1+((yrow(1+nchar)-ywas(1+nchar))/dt)**2)  ! Lorentzian Weighting
         twt=twt+wt
         do k=naux+nchar+1,ncol-1,2
            tot(k)=tot(k)+wt*(yrow(k)-ywas(k))**2/
     &      (yrow(k+1)**2+ywas(k+1)**2)
         end do
         do k=1+nchar,ncol
            ywas(k)=yrow(k)
         end do
      end do         ! do irow=1,9999999
      stop ' irow exceeded 9999999'
99    close(lunr)

c  Output results to the screen.
      do k=naux+nchar+1,ncol-1,2
         write(*,'(i4,2x,a,f9.3)') k,headarr(k),sqrt(tot(k)/twt)
      end do

      stop
      end
