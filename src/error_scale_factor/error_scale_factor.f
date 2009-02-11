c  Program: error_scale_factor.f
c
c  Purpose: To determine the factor by which the error bars
c  returned by GFIT exceed the true precision of the measurements.
c  in the selected runlog.vav.ada file
c
c
c  Input Files:
c       runlog.vav.ada 
c
c  Output Files:
c       esf_runlog.vav.ada.out
c       
c  Due to errors/ommissions in the spectroscopy and perhaps also
c  in the modeling of the instrument response (e.g. continuum, ILS),
c  the spectral fitting residuals produced by GFIT tend to be dominated
c  by systematic artifacts, not random noise. These artifacts are generally
c  the same from spectrum to spectrum, so an error calculated from the
c  rms spectral fit will generally overestimate the true precision.
c
c  This program computes the difference in the column amounts between
c  sucessive spectra and compares this with the RSS errors. The ratio
c  should average out to 1.414 if the error bars are consistent with
c  the scatter of the measured values. The weighted average value of
c  this ratio is calculated for each gas, the weights being given by
c  a Lorentzian shaped function with a half width of 5 minutes.
c
      implicit none
      integer*4 lunr,luns,ncoml,ncol,mcol,kcol,j,
     & lnbc,irow,naux,nrow,li,k
      parameter (lunr=14,luns=15,mcol=150)
      character header*800,headarr(mcol)*20,gggdir*80,
     & inputfile*40,version*60
      real*8 yrow(mcol),ywas(mcol),dt,wt,twt,tot(mcol)

      version=
     & ' apply_insitu_correction    Version 1.0.3   2009-02-01   GCT'

      call getenv('GGGPATH',gggdir)

c      open(luns,file=gggdir(:lnbc(gggdir))//'/tccon/corrections.dat',
c     & status='old')
c      read(luns,*)ncoml,ncol
c      do k=2,ncoml
c         read(luns,*)
c      end do
c      do k=1,mgas
c         read(luns,*,end=88) gasname(k),adcf(k),aicf(k)
c      end do
c      stop 'increase parameter MGAS'
c88    ngas=k-1

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
      if(kcol.ne.ncol ) stop 'ncol/kcol mismatch'

c  Read each day of data into memory and divide XGas values by the
c  appropriate correction factors.
      dt=0.00001  ! years (5 minutes)
      twt=0.0
      do k=1,ncol
         tot(k)=0.0
      end do
      read(lunr,*) (ywas(j),j=1,ncol)
      do irow=1,9999999
         read(lunr,*,end=99) (yrow(j),j=1,ncol)
         wt=1/(1+((yrow(1)-ywas(1))/dt)**2)
         do k=naux+1,ncol-1,2
            twt=twt+wt
            tot(k)=tot(k)+wt*(yrow(k)-ywas(k))**2/
     &      (yrow(k+1)**2+ywas(k+1)**2)
         end do
         do k=1,ncol
            ywas(k)=yrow(k)
         end do
      end do         ! do irow=1,9999999
      stop ' irow exceeded 9999999'
99    close(lunr)

      do k=naux+1,ncol-1,2
         write(*,*)k,headarr(k),sqrt(tot(k)/twt/2)
      end do

      stop
      end
