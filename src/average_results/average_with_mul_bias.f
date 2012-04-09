      subroutine average_with_mul_bias
     & (ymiss,nrow,ncol,yobs,yerr,
     & ybar,eybar,scal,escal,rew,cew,tew)

c  Performs a weighted average of a matrix of measurements
c  after removing any column-dependent multiplicate bias.
c
c  Decomposes the Nrow x Ncol matrix of data values into the product
c  of a row average (YBAR) and a column-dependent scaling (SCAL).
c     yobs(i,j) ~ ybar(i) * scal(j)
c
c  For example, yobs(i,j) might represent a column abundance
c  retrieved from fitting the i'th spectrum in the j'th window.
c  It is reasonable to assume that:
c    - different windows may produce systematically different
c      columns due to spectroscopic line intensity errors.
c    - different spectra have slightly different column amounts
c      due to atmospheric drifts/changes
c  in which case the entire Nrow x Ncol matrix of column amounts
c  can be boiled down to a Nrow-vector of average column amounts
c  and a Ncol-vector of window-dependent scale factors.
c  This not only represents a large reduction in data to be
c  stored, from Nrow*Ncol to Nrow+Ncol, but the resulting vectors
c  are far more meaningful, useful and less noisy.
c
c  Inputs:
c      ymiss            R*4  Missing value
c      nrow             I*4  1st dimension of YOBS (number of rows/spectra)
c      ncol             I*4  2nd dimension of YOBS (number of columns/windows)
c      yobs(nrow,ncol)  R*4  the data values
c      yerr(nrow,ncol)  R*4  the data value uncertainties
c
c  Outputs:
c      ybar(nrow)       R*4  average values of the row/spectrum
c      eybar(nrow)      R*4  average uncertainties for the row/spectrum
c      scal(ncol)       R*4  column/window scale factor
c      escal(ncol)      R*4  column/window scaling uncertainty
c      rew(nrow)        R*4  Row Error Weight (SQRT(CHI2(irow)/Ncol))
c      cew(ncol)        R*4  Scale Error Weight (SQRT(CHI2(jcol)/Nrow))
c      tew              R*4  Error Weight (SQRT(CHI2/Ncol/Nrow))
c
c  Minimizes:
c   CHI2 = Sumi Sumj [(yobs(i,j)-ybar(i)*scal(j))/yerr(i,j)]^2 + ww*Sumj(scal(j)-1)^2
c
c  Solution is obtained by differentiating the above equation
c  with respect to each element of ybar(i) and scal(j) yielding:
c   Sumj scal(j)*(yobs(i,j)-ybar(i)*scal(j))/yerr(i,j)^2 = 0  i=1,nrow
c   Sumi ybar(i)*(yobs(i,j)-ybar(i)*scal(j))/yerr(i,j)^2 + ww*(1-scal(j) = 0  j=1,ncol
c
c  which can be reorganized as
c       Sumj scal(j)*yobs(i,j)/yerr(i,j)^2  =      Sumj ybar(i)*scal(j)^2/yerr(i,j)^2  i=1,nrow
c   [ww+Sumi ybar(i)*yobs(i,j)/yerr(i,j)^2] = scal(j)*[ww+Sumi ybar(i)^2)/yerr(i,j)^2]  j=1,ncol
c
c  which can be simplified to
c    ybar(i) =    Sumj scal(j)*yobs(i,j)/yerr(i,j)^2 /    Sumj scal(j)^2/yerr(i,j)^2  i=1,nrow
c    scal(j) = ww+Sumi ybar(i)*yobs(i,j)/yerr(i,j)^2 / ww+Sumi ybar(i)^2/yerr(i,j)^2  j=1,ncol
c
c  Note that ybar(i) depends on scal(j) and vice versa. So iteration is generally needed. 
c
c  The total error weight (TEW) is the square root of:
c  the CHI2 value divided by the total number of points.
c  If the error estimates (YERR) are consistent with the
c  scatter of the data values (YOBS), TEW should be ~1.
c  If TEW>1, the scatter of the data values exceeds the error bars.
c  REW(Irow) is the same thing, but for a particular row.
c  CEW(JCOL) is the same thing, but for a particular column.
c    
c  Note that the standard errors EYBAR and ESCAL depend only on the
c  original error estimates, YERR, and not at all on the goodness
c  of fit (assuming that the SCAL values are fairly close to 1).
c  The factors REW(irow), CEW(jcol), and TEW tell you by what factor
c  the scatter of the observation (YOBS) about the mean (YBAR*SCAL)
c  exceed the errors YERR, on average.
c  There are sound arguments for multiplying EYBAR and ESCAL by TEW.
c
c  Known Problems:
c  1)  For large values of ww, iteration converges very slowly.

      implicit none
      include "params.f"

      integer*4 nrow,irow,ncol,jcol,jit,
     & nval_irow,nval_jcol,nval
      real*4 ymiss,ww,
     & yobs(nrow,ncol),yerr(nrow,ncol),
     & ybar(nrow),eybar(nrow),
     & scal(ncol),escal(ncol),rew(nrow),cew(ncol),tew,twas
      real*8 num_irow,den_irow,num_jcol,den_jcol,wt,
     & chi2_irow,chi2_jcol,chi2
      parameter (ww=0.04)   !  inverse a priori uncertainty (variance) on SCAL = 1 +/-5 

c  Check for illegal NCOL/NROW values.
      if(ncol.le.0) stop 'average_with_mul_bias:  NCOL <=0'
      if(nrow.le.0) stop 'average_with_mul_bias:  NROW <=0'

c  Initialize SCAL to a priori value.
      do jcol=1,ncol
         scal(jcol)=1.0
         escal(jcol)=1/sqrt(ww)
      end do

c  Save some time if NCOL=1
      if(ncol.eq.1) then
         call vmov(yobs,1,ybar,1,nrow)
         call vmov(yerr,1,eybar,1,nrow)
c         write(*,*)'ncol=1  Copying YOBS to YBAR'
         return
      endif

c  Now do the non-trivial cases
      twas=1.E+36
      do jit=0,mit

c  Determine YBAR, EYBAR, and REW
      do irow=1,nrow
         num_irow=0.0
         den_irow=0.0
         chi2_irow=0.0d0
         nval_irow=0.0d0
         do jcol=1,ncol
            if(yerr(irow,jcol).ne.ymiss.and.yerr(irow,jcol).ne.0.) then
               wt=scal(jcol)/yerr(irow,jcol)
               num_irow=num_irow+wt*yobs(irow,jcol)/yerr(irow,jcol)
               den_irow=den_irow+wt**2
               nval_irow=nval_irow+1
               chi2_irow=chi2_irow+
     &      ((yobs(irow,jcol)-ybar(irow)*scal(jcol))/yerr(irow,jcol))**2
            endif
         end do
         if (nval_irow.gt.0) then
            rew(irow)=sqrt(chi2_irow/nval_irow)
         else
            rew(irow)=ymiss
         endif
         if(den_irow.eq.0.0) then
            ybar(irow)=ymiss
            eybar(irow)=ymiss
         else
            ybar(irow)=num_irow/den_irow
            eybar(irow)=1/sqrt(den_irow)
         endif
      end do

c  Determine SCAL, ESCAL, and CEW
      chi2=0.0
      nval=0
      do jcol=1,ncol
         num_jcol=ww
         den_jcol=ww 
         chi2_jcol=ww*(1-scal(jcol))**2
         nval_jcol=0
         do irow=1,nrow
            if(yerr(irow,jcol).ne.ymiss.and.yerr(irow,jcol).ne.0.) then
               wt=ybar(irow)/yerr(irow,jcol)
               num_jcol=num_jcol+wt*yobs(irow,jcol)/yerr(irow,jcol)
               den_jcol=den_jcol+wt**2
               nval_jcol=nval_jcol+1
               chi2_jcol=chi2_jcol+
     &      ((yobs(irow,jcol)-ybar(irow)*scal(jcol))/yerr(irow,jcol))**2
            endif
         end do
         if(nval_jcol.gt.0) then
            cew(jcol)=sqrt(chi2_jcol/nval_jcol)
         else
            cew(jcol)=ymiss
         endif
         nval=nval+nval_jcol
         chi2=chi2+chi2_jcol
         if(den_jcol.eq.0.0) then
            write(*,*)' No spectra for this window'
            scal(jcol)=ymiss
            escal(jcol)=ymiss
         else 
            scal(jcol)=num_jcol/den_jcol
            escal(jcol)=1/sqrt(den_jcol)
         endif
      end do
      tew=sqrt(chi2/nval)
c      write(*,*)jit,twas,tew
      if( tew .ge. twas ) exit  ! fit failed to improve
      if(jit.eq.mit) write(*,*)
     & 'average_with_mul_bias failed to converge',twas,tew
      twas=tew
      end do  !  jit=1,mit
c
c  Scale EYBAR by MAX(1,TEW)
c  Not scaling EYBAR would imply that the GFIT-supplied error bars are correct
c  Scaling EYBAR by TEW would imply that the residuals represent the measurement incertainties
      call vmul(eybar,1,amax1(1.0,tew),0,eybar,1,nrow)
      return
      end
