      subroutine average_with_mul_bias
     & (mit,ymiss,nrow,ncol,yobs,yerr,
     & ybar,eybar,scal,serr,rew,cew,tew)

c  Performs a weighted average for each row of a matrix of measurements
c  after removing any column-dependent multiplicate biases.
c
c  Inputs:
c     mit             I*4  Max # of iterations
c     ymiss           R*4  Missing value
c     nrow            I*4  1st dimension of YOBS (number of rows/spectra)
c     ncol            I*4  2nd dimension of YOBS (number of columns/windows)
c     yobs(nrow,ncol) R*4  the data values
c     yerr(nrow,ncol) R*4  the data value uncertainties
c     scal(ncol)      R*4  A Priori scale factor values
c     serr(ncol)      R*4  A Priori scale factor uncertainties
c
c  Outputs:
c     ybar(nrow)      R*4  average values of the row/spectrum
c     eybar(nrow)     R*4  average uncertainties for the row/spectrum
c     scal(ncol)      R*4  New column/window scale factor
c     serr(ncol)      R*4  New column/window scaling uncertainty
c     rew(nrow)       R*4  Row Error Weight (SQRT(CHI2(irow)/Ncol))
c     cew(ncol)       R*4  Column Error Weight (SQRT(CHI2(jcol)/Nrow))
c     tew             R*4  Total Error Weight (SQRT(CHI2/Ncol/Nrow))
c
c  Decomposes a Nrow x Ncol matrix of data values into the product
c  of a row average (YBAR) and a column-dependent scaling (SCAL).
c     yobs(i,j) ~ ybar(i) * scal(j)
c
c  For example, yobs(i,j) might represent a column abundance
c  retrieved from fitting the i'th spectrum in the j'th window.
c  It is reasonable to assume that:
c    - different windows (columns) may produce systematically
c      different column abundances due to spectroscopic errors.
c    - different spectra (rows) have different column amounts
c      due to atmospheric drifts/changes
c  in which case the entire Nrow x Ncol matrix of column amounts
c  can be reduced to a Nrow-vector of average column amounts and
c  a Ncol-vector of window-dependent scale factors. This reduces
c  the data to be stored from Nrow*Ncol to Nrow+Ncol, with little
c  loss of information because initially there was much redundacy.
c
c  The ybar and scal vectors are determined by minimizing:
c
c  Chi^2 = Sumi Sumj [(yobs(i,j)-ybar(i)*scal(j))/yerr(i,j)]^2 +
c          Sumj[(scal(j)-aps(j))/serr(j)]^2
c
c  The upper term tries to find the ybar and scal values that best
c  fit the data. The lower term is a weak a priori constraint that
c  tries to keep the scale factors close to aps (usually ~1.0).
c
c  Solution is obtained by differentiating the above equation
c  with respect to each element of ybar(i) and scal(j) yieldingi
c  the simultaneous equations:
c
c   Sumj scal(j)*(yobs(i,j)-ybar(i)*scal(j))/yerr(i,j)^2 = 0   i=1,nrow
c
c   Sumi ybar(i)*(yobs(i,j)-ybar(i)*scal(j))/yerr(i,j)^2 + 
c    (aps(j)-scal(j))/serr(j)^2 = 0                        j=1,ncol
c
c  which can be simplified to:
c
c            Sumj[scal(j)*yobs(i,j)/yerr(i,j)^2]
c  ybar(i) = -----------------------------------         i=1,nrow
c               Sumj[scal(j)^2/yerr(i,j)^2]
c
c            aps(j)/serr(j)^2 + Sumi[ybar(i)*yobs(i,j)/yerr(i,j)^2]
c  scal(j) = ------------------------------------------------------  j=1,ncol
c               1/serr(j)^2   +    Sumi[ybar(i)^2/yerr(i,j)^2] 
c
c N.B:
c - Since ybar(i) depends on scal(j) and vice versa, iteration is needed.
c   which is another way of saying that the chi^2 equation is non-linear.
c - As serr gets smaller, or yerr bigger,  scal(j) --> aps(j)
c - Although ybar depends on scal, it is independent of serr.
c
c  The total error weight (TEW) is the square root of the CHI2 value
c  divided by the total number of points. If the error bars (YERR)
c  are consistent with the scatter of the data values (YOBS), TEW
c  should be ~1. If TEW>1, the scatter of the data values exceeds
c  the error bars.
c  REW(Irow) is the same thing, but for a particular row (spectrum).
c  CEW(JCOL) is the same thing, but for a particular column (window).
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
c  1) For small values of apserr, iteration converges very slowly,
c   resulting in warnings..

      implicit none
c      include "params.f"

      integer*4 nrow,irow,ncol,jcol,jit,mit,
     & nval_irow,nval_jcol,nval
      real*4 ymiss,sf,small,
c     & yemin,yemax,
     & yobs(nrow,ncol),yerr(nrow,ncol),
     & ybar(nrow),eybar(nrow),
     & apscal(ncol),apserr(ncol),
     & scal(ncol),serr(ncol),rew(nrow),cew(ncol),tew,twas
      real*8 num_irow,den_irow,num_jcol,den_jcol,wt,
     & chi2_irow,chi2_jcol,chi2
      logical isclose_s  ! function to compare values with floating point error

       small=1.0e-18
c      parameter (ww=0.04)   !  inverse a priori uncertainty (variance) on SCAL = 1 +/-5 

c  Check for illegal NCOL/NROW values.
      if(ncol.le.0) stop 'average_with_mul_bias:  NCOL <=0'
      if(nrow.le.0) stop 'average_with_mul_bias:  NROW <=0'

c      write(*,*)'navg=',ncol
c      do jcol=1,ncol
c         yemin=yerr(1,jcol)
c         yemax=yerr(1,jcol)
c         do irow=2,nrow
c            if(yerr(irow,jcol).lt.yemin) yemin=yerr(irow,jcol)
c            if(yerr(irow,jcol).gt.yemax) yemax=yerr(irow,jcol)
c         end do
c         write(*,*)'jcol,yemin,yemax=',jcol,yemin,yemax
c      end do

c  Set a prior SCAL values and uncertainties.
      do jcol=1,ncol
c         write(*,*) 'jcol,scal=',jcol,scal(jcol)
         apscal(jcol)=scal(jcol)
         apserr(jcol)=serr(jcol)
      end do

c  Save some time if NCOL=1
      if(ncol.eq.1) then
c         call vmov(yobs,1,ybar,1,nrow)
c         call vmov(yerr,1,eybar,1,nrow)
         do irow=1,nrow
            ybar(irow)=yobs(irow,1)
            eybar(irow)=yerr(irow,1)
         end do
c         write(*,*)'ncol=1  Copying YOBS to YBAR'
         return
      endif

c  Now do the non-trivial cases
      twas=1.E+36
      do jit=1,mit

c  Determine YBAR, EYBAR, and REW
         do irow=1,nrow
c            write(*,*) 'irow,ybar(irow)  aa=',irow,ybar(irow)
            num_irow=small
            den_irow=small
            num_irow=0.0d0
            den_irow=0.0d0
            chi2_irow=0.0d0
            nval_irow=0
            do jcol=1,ncol
c               if(yerr(irow,jcol).ne.ymiss .and.
c     &         yerr(irow,jcol).ne.0.) then
c                if( abs(yerr(irow,jcol)-ymiss).gt.small 
c     &         .and. abs(yerr(irow,jcol)).gt.small) then
               if( .not. isclose_s(yerr(irow,jcol), ymiss) ! if(yerr.ne.ymiss)
     &         .and. abs(yerr(irow,jcol)).gt.small ) then ! if(yerr.ne.0)
                  wt=scal(jcol)/yerr(irow,jcol)
                  num_irow=num_irow+wt*yobs(irow,jcol)/yerr(irow,jcol)
                  den_irow=den_irow+wt**2
                  nval_irow=nval_irow+1
                  chi2_irow=chi2_irow+
     &            ((yobs(irow,jcol)-ybar(irow)*scal(jcol))/
     &            yerr(irow,jcol))**2
               endif
            end do
            if (nval_irow.gt.0) then
               ybar(irow)=sngl(num_irow/den_irow)
               eybar(irow)=sngl(1/dsqrt(den_irow))
               rew(irow)=sngl(dsqrt(chi2_irow/nval_irow))
            else
               ybar(irow)=ymiss
               eybar(irow)=ymiss
               rew(irow)=ymiss
            endif
c            write(*,*) 'irow,ybar(irow) bb=',irow,ybar(irow),
c     &      num_irow,den_irow
         end do  ! do irow=1,nrow

c  Determine SCAL, ESCAL, and CEW
         chi2=0.0
         nval=0
         do jcol=1,ncol
            num_jcol=1/apserr(jcol)**2
            den_jcol=1/apserr(jcol)**2
c            write(*,*)'jcol,num_jcol,den_jcol 1=',jcol,num_jcol,den_jcol
            chi2_jcol=((apscal(jcol)-scal(jcol))/apserr(jcol))**2
c            write(*,*)'jit,jcol,spscal,apserr=',
c     &      jit,jcol,apscal(jcol),apserr(jcol),apscal(jcol),scal(jcol)
            nval_jcol=0
            do irow=1,nrow
c               if( .not. isclose_s(yerr(irow, jcol), ymiss)  ! if(yerr.ne.ymiss)
                if(abs(yerr(irow,jcol)-ymiss) .gt. small
     &         .and. abs(yerr(irow,jcol)) .gt. small) then  ! if(yerr.ne.0)
                  wt=ybar(irow)/yerr(irow,jcol)
                  num_jcol=num_jcol+wt*yobs(irow,jcol)/yerr(irow,jcol)
                  den_jcol=den_jcol+wt**2
                  nval_jcol=nval_jcol+1
c              write(*,*)'yy=',irow,jcol,ybar(irow),scal(jcol),
c     &        ((yobs(irow,jcol)-ybar(irow)*scal(jcol))/
c     &            yerr(irow,jcol))**2
                  chi2_jcol=chi2_jcol+
     &            ((yobs(irow,jcol)-ybar(irow)*scal(jcol))/
     &            yerr(irow,jcol))**2

                  if (((yobs(irow,jcol)-ybar(irow)*scal(jcol))/
     &            yerr(irow,jcol))**2 .gt. ymiss) then
                     stop 'chi2_jcol = Inf'
                  endif
               else
c                  write(*,*)'yerr=ymiss or 0',yerr(irow,jcol)
               endif
            end do
            if(nval_jcol.gt.0) then
               cew(jcol)=sngl(dsqrt(chi2_jcol/nval_jcol))
            else
               cew(jcol)=ymiss
            endif
c            write(*,*)'jcol,chi2_jcol=',jcol,chi2_jcol
            nval=nval+nval_jcol
            chi2=chi2+chi2_jcol
c            write(*,*)'jcol,num_jcol,den_jcol 2=',jcol,
c     &      num_jcol,den_jcol,chi2_jcol
            if(den_jcol.lt.small) then
               write(*,*)' No spectra for this window'
               scal(jcol)=ymiss
               serr(jcol)=ymiss
            else 
               scal(jcol)=sngl(num_jcol/den_jcol)
               serr(jcol)=sngl(1/dsqrt(den_jcol))
            endif
         end do          ! do jcol=1,ncol
c         write(*,*)'chi2,nval=',chi2,nval
         tew=sngl(dsqrt(chi2/nval))
         if( tew .ge. twas ) exit  ! fit failed to improve
         if(mit.gt.1 .and. jit.eq.mit) write(*,*)
     &    'average_with_mul_bias failed to converge',twas,tew
         twas=tew
      end do  !  jit=1,mit
c      write(*,*)'mit,nit,tew=',mit,jit-1,tew,chi2,nval
c
c  Scale EYBAR by MAX(1,TEW) (unless EYBAR is a fill value)
c  Not scaling EYBAR would imply that the GFIT-supplied error bars are correct
c  Scaling EYBAR by TEW would imply that the residuals represent the measurement incertainties
      sf=amax1(1.0,tew)
      do irow=1,nrow
        if(.not. isclose_s(eybar(irow), ymiss))
     &      eybar(irow)=sf*eybar(irow)
      end do


      return
      end
