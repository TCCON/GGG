      subroutine wlsfit(mrow,nmp,yydu,nfp,bfdu,apx,apu,rnorm)
c  Performs a weighted least-squares fit of a linear multi-variate
c  function f(i) to the data yy(i) +/- uy(i)  by minimizing
c     Chi2 = Sum_i [((yy(i)-f(i))/uy(i)]**2
c  where y(i) are the measurements having uncertainties uy(i)
c  and f(i) is a linear model to be fitted to the data.
c     f(i) = Sum_j[CX(j)*BF(i,j)]
c  where CX(j) is the unknown parameters to be retrieved
c  and BF(i,j) are the basis functions.

c Note that the YY vector and the BF matrix passed to this subroutine
c must have already be divided by the measurement uncertainties UY
c
c  Inputs:
c       NPTS           I*4  Number of measured points
c       YYDU(NMP)      R*4  Measured data values divided by their uncertainties
c       NFP            I*4  Number of fitted parameters
c       BFDU(NMP,NFP)  R*4  Array of basis functions divided by the measurement uncertainties
c       APX(NFP)       R*4  A priori estimate of solution vector
c       APU(NFP)       R*4  Uncertainty of a priori estimate
c
c Outputs:
c       (YYDU(k),k=1,nfp)          the retrieved coefficients (CX)
c       (sqrt(BFDU(k,k)),k=1,nfp)  the uncertainties in CX 
c       RNORM                      R*4  Predicted L2-Norm (chi-squared = RNORM**2)
c
c   RNORM/SQRT(NMP)  is the factor by which the scatter of the points exceeds their error bars

      implicit none
      integer*4 j,k,mrow,nmp,mfp,nfp,krank,ierr
      parameter (mfp=20)
      integer*4 ip(mfp)
      real*4 yydu(mrow),bfdu(mrow,nfp),apx(nfp),apu(nfp),wk(mfp),
     & tau,rnorm,var

      tau=1.0E-06
      if(nfp.gt.mfp) stop ' LSFIT: NFP > MFP'
      if(nmp+nfp.gt.mrow) stop ' LSFIT: NMP+NFP > MROW'

c  Add a priori info to bottom of BFDU and YY arrays
      do k=1,nfp
         do j=1,nfp
            bfdu(nmp+k,j)=0.0
         end do
         bfdu(nmp+k,k)=1/apu(k)
         yydu(nmp+k)=apx(k)/apu(k)
      end do

c  Perform QR decomposition of matrix BFDU and vector YY
      call shfti(bfdu,mrow,nmp+nfp,nfp,yydu,nmp+nfp,
     & 1,tau,krank,rnorm,wk,ip)
      if(krank.lt.nfp) write(*,*) 'Warning: rank-deficiency:',krank,nfp
      var=1.0  ! Since we've already divided everything by the measurement uncertainties.
      call scov2(bfdu,mrow,nfp,ip,var,ierr)
      if(ierr.ne.0) write(*,*)'ierr=',ierr
      return
      end 
