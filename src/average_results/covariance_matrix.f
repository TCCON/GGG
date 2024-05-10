      subroutine covariance_matrix
     & (ymiss,nexp,ncol,nrow,yobs,yerr,covmat)

c  Computes a matrix of weighted covariances
c  between the data in the various columns of YOBS.
c  covmat(k,j) = SUM_i wt.[yobs(i,k)-ybar(k)].[yobs(i,j)-ybar(j)]
c
c  where yobs(i,k,i) might represent a column abundance
c  retrieved from fitting the i'th spectrum in the k'th window.
c  and yerr(k,i) represents its uncertainty.
c  ybar(k) is the weighted average column abundance
c
c  wt = [1.0d0/yerr(i,k)/yerr(i,k)]^NEXP
c
c  Points with large YERR can be de-weighted by choosing NEXP>0
c  In this case the returned matrix is not the true covariances.
c  Using NEXP=0 disregards uncertainties and should be used only
c  for testing purposes, for example validating against other
c  calculations of the correlation coefficients (e.g. MATLAB).
c
c  Inputs:
c      ymiss            R*4  Missing value
c      nexp             I*4  Exponent of weights
c      ncol             I*4  1st dimension of YOBS (number of columns/windows)
c      nrow             I*4  2nd dimension of YOBS (number of rows/spectra)
c      yobs(nrow,ncol)  R*4  data values
c      yerr(nrow,ncol)  R*4  data value uncertainties
c
c  Outputs:
c      covmat(ncol,ncol)   R*8  weighted covariance matrix
c
      implicit none
      integer*4 nrow,irow,ncol,kcol,jcol,nexp
      real*4 ymiss,small,yobs(nrow,ncol),yerr(nrow,ncol),
     & ymin,ymax,emin,emax,
     & ybar(ncol)
      real*8 wt,toty,totw,totwc,covmat(ncol,ncol)

c  Check for illegal NCOL/NROW values.
      if(ncol.le.0) stop 'weighted_covariance_matrix:  NCOL <=0'
      if(nrow.le.0) stop 'weighted_covariance_matrix:  NROW <=0'

c  Compute weighted mean value of each column (window)
      small=1.0E-18
      ymin=1.E+36
      emin=1.E+36
      ymax=-1.E+36
      emax=-1.E+36
      do kcol=1,ncol
         totw=0.0d0
         toty=0.0d0
         do irow=1,nrow
            if(abs(yobs(irow,kcol)-ymiss).gt.small) then
               wt=dabs(1.0d0/yerr(irow,kcol))**nexp
               totw=totw+wt
               toty=toty+wt*yobs(irow,kcol)
               if(yobs(irow,kcol).gt.ymax) ymax=yobs(irow,kcol)
               if(yobs(irow,kcol).lt.ymin) ymin=yobs(irow,kcol)
               if(yerr(irow,kcol).gt.emax) emax=yerr(irow,kcol)
               if(yerr(irow,kcol).lt.emin) emin=yerr(irow,kcol)
            endif
         end do
         ybar(kcol)=sngl(toty/totw)
c         write(*,*) 'kcol, ybar = ',kcol, ybar(kcol),toty,totw
      end do
c      write(*,*)'ymin, emin=',ymin,emin
c      write(*,*)'ymax, emax=',ymax,emax

c  Computed weighted covariances
      do kcol=1,ncol
         do jcol=1,kcol
            totwc=0.0d0
            do irow=1,nrow
c               if(yobs(irow,kcol).ne.ymiss) then
c               if(yobs(irow,jcol).ne.ymiss) then
               wt=dabs(1.0d0/yerr(irow,kcol)/yerr(irow,jcol))**nexp
               totwc=totwc+wt*
     &         dble(yobs(irow,kcol)-ybar(kcol))*
     &             (yobs(irow,jcol)-ybar(jcol))
c               endif
c               endif
            end do  ! irow=1,nrow
            covmat(kcol,jcol)=totwc
            covmat(jcol,kcol)=totwc
c            write(*,*)kcol,jcol,totwc
         end do   ! jcol=1,ncol
      end do   ! kcol=1,ncol
      
      return
      end
