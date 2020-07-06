      subroutine average_with_add_bias
     & (mit,ymiss,nrow,ncol,yobs,yerr,
     & ybar,eybar,bias,ebias,rew,cew,tew)
c
c  Performs a weighted average of a series of measurements
c  after removing any additive bias.
c
c  Decomposes a Nrow x Ncol matrix of data values into the sum of a
c  row average (YBAR) and a column-dependent bias or offset (BIAS).
c     yobs(i,j) ~ ybar(i) + bias(j)
c
c  For example, yobs(i,j) might represent a frequency stretch
c  retrieved from fitting the i'th spectrum in the j'th window.
c  It is reasonable to assume that:
c    - different windows may produce systematically different
c      stretches due to spectroscopic line position errors.
c    - different spectra have slightly different frequency
c      calibrations due to instrumental drifts/changes
c  in which case the entire Nrow x Ncol matrix of frequency shifts
c  can be boiled down to a Nrow-vector of spectral frequency offsets
c  and a Ncol-vector of window-dependent frequency offsets/biases.
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
c      ybar(nrow)       R*4  row-average values
c      eybar(nrow)      R*4  row-average uncertainties
c      bias(ncol)       R*4  column bias
c      ebias(ncol)      R*4  column bias uncertainty
c      rew(nrow)        R*4  Column Error Weight = SQRT(CHI2(jcol)/Nrow)
c      cew(ncol)        R*4  Column Error Weight = SQRT(CHI2(jcol)/Nrow)
c      tew              R*4  Total Error Weight  = SQRT(CHI2/(Nrow*Ncol))
c
c  Minimizes
c     CHI2 = Sumi Sumj [(yobs(i,j)-ybar(i)-bias(j))/yerr(i,j)]^2 +ww*Sumj(bias(j)^2)
c
c  Solution is obtained by differentiating the above equation
c  with respect to each element of ybar(i) and bias(j) yielding:
c     Sumj (yobs(i,j)-ybar(i)-bias(j))/yerr(i,j)^2 = 0  i=1,nrow
c     Sumi (yobs(i,j)-ybar(i)-bias(j))/yerr(i,j)^2 - ww*bias(j) = 0  j=1,ncol
c
c  which can be reorganized as
c     ybar(i) = Sumj (yobs(i,j)-bias(j))/yerr(i,j)^2 / Sumj 1/yerr(i,j)^2  i=1,nrow
c     bias(j) = Sumi (yobs(i,j)-ybar(i))/yerr(i,j)^2 / ww+Sumi 1/yerr(i,j)^2  j=1,ncol
c
c  Note that ybar(i) depends on bias(j) and vice versa.  So that iteration is generally needed.
c  When the uncertainties yerr(i,j) are constant, however, I have found empirically that
c  iteration is not necessary (for reasons that I don't fully understand).
c
c  The total error weight (TEW) is the square root of:
c    the CHI2 value divided by the total number of points.
c  If the error estimates (YERR) are consistent with the
c  scatter of the data values (YOBS), TEW should be ~1.
c  If TEW>1, the scatter of the data values exceeds the error bars.
c  CEW(JCOL), is the same thing, but for a particular column.
c
c  Note that the standard errors EYBAR and EBIAS depend only on the
c  original error estimates, YERR, and not at all on the goodness
c  of fit (assuming that the BIAS values are fairly close to 1).
c  The factors REW(irow), CEW(jcol), and TEW tell you by what factor
c  the scatter of the observation (YOBS) about the mean (YBAR*BIAS)
c  exceed the errors YERR, on average.
c  There are sound arguments for multiplying EYBAR and EBIAS by TEW.
c
c  Known Problems:
c  1)  For large values of ww, iteration converges very slowly.

      implicit none
c      include "params.f"

      integer*4 nrow,irow,ncol,jcol,jit,mit,
     & nval_irow,nval_jcol,nval
      real*4 ymiss,ww,small,
     & yobs(nrow,ncol),yerr(nrow,ncol),
     & ybar(nrow),eybar(nrow),
     & bias(ncol),ebias(ncol),rew(nrow),cew(ncol),tew,twas
      real*8 num_irow,den_irow,num_jcol,den_jcol,wt,
     & chi2_irow,chi2_jcol,chi2

      small=1.0E-18
c  Check for illegal NCOL/NROW values.
      if(ncol.le.0) stop 'average_with_add_bias:  NCOL <=0'
      if(nrow.le.0) stop 'average_with_add_bias:  NROW <=0'

c  Initialize BIAS to a priori value.
      do jcol=1,ncol
         bias(jcol)=0.0
         ebias(jcol)=0.0
      end do

c  Save some time if NCOL=1
      if(ncol.eq.1) then
         call vmov(yobs,1,ybar,1,nrow)
         call vmov(yerr,1,eybar,1,nrow)
c         write(*,*)'ncol=1  Copying YOBS to YBAR'
         return
      endif

c  Now for the non-trivial cases.
      ww=0.0E-16  ! inverse a priori variance of bias
      twas=1.E+36
      do jit=1,mit

c  Determine YBAR, EYBAR, and REW
         do irow=1,nrow
            num_irow=0.0d0
            den_irow=0.0d0
            chi2_irow=0.0d0
            nval_irow=0
            do jcol=1,ncol
               if(abs(yerr(irow,jcol)-ymiss).gt.small .and.
     &         abs(yerr(irow,jcol)).gt.small) then
                  wt=1.0d0/yerr(irow,jcol)**2
                  num_irow=num_irow+(yobs(irow,jcol)-bias(jcol))*wt
                  den_irow=den_irow+wt
                  nval_irow=nval_irow+1
                  chi2_irow=chi2_irow+
     &            wt*(yobs(irow,jcol)-ybar(irow)-bias(jcol))**2
               endif
            end do
            rew(irow)=sngl(dsqrt(chi2_irow/nval_irow))
            if(abs(den_irow).lt.small) then
               ybar(irow)=ymiss
               eybar(irow)=ymiss
            else
               ybar(irow)=sngl(num_irow/den_irow)
               eybar(irow)=sngl(1/dsqrt(den_irow))
            endif
         end do

c  Determine BIAS, EBIAS, and CEW
         chi2=0.0
         nval=0
         do jcol=1,ncol
            num_jcol=0.0
            den_jcol=ww      ! a priori bias = 0 +/- 1/SQRT(WW)
            chi2_jcol=ww*bias(jcol)**2
            nval_jcol=0
            do irow=1,nrow
               if(abs(yerr(irow,jcol)-ymiss).gt.small .and.
     &         abs(yerr(irow,jcol)).gt.small) then
                  wt=1.0d0/yerr(irow,jcol)**2
                  num_jcol=num_jcol+(yobs(irow,jcol)-ybar(irow))*wt
                  den_jcol=den_jcol+wt
                  nval_jcol=nval_jcol+1
                  chi2_jcol=chi2_jcol+
     &            wt*(yobs(irow,jcol)-ybar(irow)-bias(jcol))**2
               endif
            end do
            cew(jcol)=sngl(dsqrt(chi2_jcol/nval_jcol))
            nval=nval+nval_jcol
            chi2=chi2+chi2_jcol
            if(abs(den_jcol).lt.small) then
               write(*,*)' No spectra for this window'
               bias(jcol)=ymiss
               ebias(jcol)=ymiss
            else 
               bias(jcol)=sngl(num_jcol/den_jcol)
               ebias(jcol)=sngl(1/dsqrt(den_jcol))
            endif
         end do
         tew=sngl(dsqrt(chi2/nval))
         if( tew .ge. twas ) exit  ! fit failed to improve
         twas=tew
      end do  !  jit=1,mit
      if(jit.gt.mit) write(*,*)
     & 'average_with_mul_bias failed to converge',twas,tew
c
c  Scale EYBAR by MAX(1,TEW)
c  Not scaling EYBAR would imply that the GFIT-supplied error bars are correct
c  Scaling EYBAR by TEW would imply that the residuals represent the measurement incertainties
      call vmul(eybar,1,amax1(1.0,tew),0,eybar,1,nrow)
      return
      end

