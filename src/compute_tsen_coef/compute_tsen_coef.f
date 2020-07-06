c  Program compute_tsen_coef.f
c
c  Parameterizes the T-sensitivity of all the TCCON windows
c  by fitting a function of the form
c    Yobs(l) = 1+c1+c2.Sqrt[airmass(l)]+c3.h2o(l)+c4.Tsurf(l)
c  to the inputted column ratios, from synthetic spectra
c  computed with and without the 7/5/3/1 K T-perturbation.
c  l is an index over spectra/rows of the .vsw file
c  Tpert is the temperature perturbation applied to the .mod
c  file in generating the fitted spectra.
c  Computes T-correction regression coefficients (c1,c2,c3,c4) by
c  solving
c     A.X=Yobs
c  where
c  A is the (MxN) matrix of basis functions (airmass, H2O, T)
c  X is the (NxK) matrix of unknown regression coefficients
c  Yobs is the (MxK) table of perturbed VSF values
c  M (NROW) is the number of spectra (measured VSF values per window)
c  N (NFC)  is the number of fitted coefficients
c  K (NWIN) is the number of windows (right-hand sides)
c
c Program is run by typing:
c       ./compute_tcorr_coeff tsen_ratio_vsw.out N
c  where the .out file contains Yobs, the T-perturbed VSF values
c  for each window and spectrum. The "N" represents NFC, the number
c  of fitted coefficients.
c
c  Input files:
c    tsen_ratio.out  ! File containing VSFs resulting from T-perturbation
c
c Output files:  all written to local directory
c    tsen.nam0    ! List of window names and RMS for the do-nothing case
c    tsen.rmsN    ! RMS of fit to VSFs for each window (for N=NFC)
c    tsen.coefN   ! Regression coefficients (for N=NFC)
c
c  tsen.coeff files must be copied to $GGGPATH/tsen_coefs/ 
c  before they can be used by apply_tcor.

      implicit none

      integer irow,i,j,k,mlabel,ipasza,
     & lunr_xsw,lunw_coef,lunw_nam,lunw_rms,
     & mfc,nfc,krank,
     & mwin,nwin,mauxcol,
     & mrow,
     & nss,
     & nauxcol,nlhead,ncol,nrow

      parameter (lunr_xsw=51)      ! input file (.xsw)
      parameter (lunw_coef=52)     ! output file (.coef)
      parameter (lunw_nam=53)      ! window names output file (.rms)
      parameter (lunw_rms=54)      ! output file (.rms)
      parameter (mauxcol=30)       !
      parameter (mwin=240)         ! Total number of columns/windows
      parameter (mrow=50)          ! Max number of output records/spectra
      parameter (mlabel=18000)     ! Max Number of column lable characters
      parameter (mfc=5)            ! Number of fitted coefficients

      integer ip(mfc)

      character
     & version*64,
     & tsen_ratio_version*64,
     & gfit_version*64,
     & gsetup_version*64,
     & collate_version*64,
     & tsen_ratio_file*80,
     & collabel*(mlabel),
     & vs(2*mwin+mauxcol)*10,
     & col1(mrow)*1,
     & cdum*1,
     & spectrum*57, 
     & data_fmt*48

      real*8 year(mrow)

      real*4
     & a(mrow,mfc),work(mfc),rnorm(mwin),tau,
     & airmass,wvcol,tout,d2r,
     & ymiss,
     & yaux(mrow,mauxcol),
     & yobs(mrow,mwin),yerr(mrow,mwin)

      version=
     & ' compute_tsen_coef.f      Version 1.05        2019-07-12   GCT'

      d2r=3.14159265/180.0
      tau=1.E-5

c  Read input info from command line
      if (iargc() == 2) then
         call getarg(1, tsen_ratio_file)
         call getarg(2, cdum)
         read(cdum,*) nfc
         if(nfc.gt.mfc) stop 'compute_tsen_coef: NFC > MFC'
      else
         stop 'Usage: compute_tsen_coef tsen_ratio.out nfc'
      endif

c  Read the entire contents of the .xsw disk file
      open(lunr_xsw,file=tsen_ratio_file,status='old')
      read(lunr_xsw,*) nlhead,ncol,nrow,nauxcol
      if(nlhead.ne.8) stop 'compute_tsen_coef: nlhead .ne. 8'
      if(nrow.gt.mrow) stop 'compute_tsen_coef: increase parameter mrow'
      read(lunr_xsw,'(a)') tsen_ratio_version
      read(lunr_xsw,'(a)') collate_version
      read(lunr_xsw,'(a)') gfit_version
      read(lunr_xsw,'(a)') gsetup_version
      read(lunr_xsw,'(8x,e12.5)') ymiss
      read(lunr_xsw,'(7x,a)') data_fmt
      read(lunr_xsw,'(a)') collabel
      call substr(collabel,vs,2*mwin+mauxcol,nss)
      if(nss.ne.ncol) then
         write(*,*)' nss, ncol = ', nss,ncol
         write(*,*)' nss .ne. ncol '
      endif
      ipasza=0
      do i=1,nauxcol
         if(vs(i).eq.'asza') ipasza=i
      end do

      nwin=(ncol-nauxcol)/2
      write(*,*) 'nwin, mwin',nwin,mwin
      if(nwin.gt.mwin) then
         write(*,*) 'nwin, mwin',nwin,mwin
         stop 'compute_tsen_coef: nwin > mwin'
      endif
      write(*,*) 'data fmt =  ',data_fmt

c  Open output file
      open(lunw_nam,file='tsen_window.nam',status='unknown')
      do k=1,nwin
         write(lunw_nam,'(a10)') vs(nauxcol+2*k-1)
      end do
      close(lunw_nam)

      do irow=1,mrow
         if (index(collabel, 'spectrum') .gt. 0) then
            read(lunr_xsw,data_fmt,end=99) spectrum,
     &      col1(irow),year(irow),(yaux(irow,k),k=3,nauxcol),
     $      (yobs(irow,k),yerr(irow,k),k=1,nwin)
         else
            read(lunr_xsw,data_fmt,end=99)
     $      col1(irow),year(irow),(yaux(irow,k),k=2,nauxcol),
     $      (yobs(irow,k),yerr(irow,k),k=1,nwin)
         endif

c Subtract 1.0 from the ratios
         do k=1,nwin
            yobs(irow,k)=yobs(irow,k)-1.0
         end do

c  Compute Airmass basis functions
         airmass=1.01/(0.02+cos(d2r*yaux(irow,ipasza)))

c  Temperatures and H2O columns used in simulations
         if(irow.gt.30) then
            tout=305.-273.15
            wvcol=1E+23
         elseif(irow.gt.25) then
            tout=305.-273.15
            wvcol=3E+22
         elseif(irow.gt.20) then
            tout=285.-273.15
            wvcol=1E+23
         elseif(irow.gt.15) then
            tout=285.-273.15
            wvcol=3E+22
         elseif(irow.gt.10) then
            tout=285.-273.15
            wvcol=1E+22
         elseif(irow.gt.5) then
            tout=265.-273.15
            wvcol=3E+22
         else
            tout=265.-273.15
            wvcol=1E+22
         endif
         if(nfc.ge.1) a(irow,1)=1.
         if(nfc.ge.2) a(irow,2)=sqrt(airmass)
         if(nfc.ge.3) a(irow,3)=wvcol*1.E-22
         if(nfc.ge.4) a(irow,4)=tout
         if(nfc.ge.5) a(irow,5)=wvcol*1.E-22*airmass
         write(*,*) irow,(a(irow,j),j=1,nfc),1+yobs(irow,3)

      end do  !  irow=1,mrow
99    close(lunr_xsw)
      if(nrow.ne.irow-1) then
         write(*,*) 'nrow, irow-1=',nrow,irow-1
         stop 'compute_tsen_coef: NROW mismatch'
      endif
c
c  Do the linear regression, with multiple (NWIN) RH sides.
      call SHFTI (a,mrow,nrow,nfc,yobs,mrow,nwin,tau,KRANK,
     & rnorm,work,ip)

c  Yobs now contains X, the regression coefficients
      open(lunw_rms, file='tsen.rms'//cdum,status='unknown')
      open(lunw_coef,file='tsen.coef'//cdum,status='unknown')
      write(lunw_coef,*) 7,nfc+1
      write(lunw_coef,'(a)') version
      write(lunw_coef,'(a)') tsen_ratio_version
      write(lunw_coef,'(a)') collate_version
      write(lunw_coef,'(a)') gfit_version
      write(lunw_coef,'(a)') gsetup_version
      write(lunw_coef,'(a6,6(a9,i1))') 'window',('        c',i,i=1,nfc)
      do k=1,nwin
         write(lunw_coef,'(a10,6(1x,f9.6))') vs(nauxcol+2*k-1),
     &   (yobs(j,k),j=1,nfc)
         write(lunw_rms,'(f7.5)') rnorm(k)/sqrt(float(nrow))
      end do
      close(lunw_coef)
      close(lunw_rms)

      stop
      end
