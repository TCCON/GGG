      function log1pxox(x)

c  Computes loge(1+x)/x without risk of underflow or zero-division.
c  For x > 0.1, loge(1+x)/x is computed directly.
c  For x < 0.1, LOG1PXOX is approximated by the series:
c  1 - x/2 + x^2/3 - x^3/4 + x^4/5 - x^5/6 + ...
c  For single-precision work,  6 terms should be taken.
c  For double-precision work, 13 terms should be taken.
c
c    J       log1pxox
c  --------------------------------------------------------
c    6     1/6
c    5     1/5-x/6
c    4     1/4-x(1/5-x/6)
c    3     1/3-x(1/4-x(1/5-x/6))
c    2     1/2-x(1/3-x(1/4-x(1/5-x/6)))
c    1     1/1-x(1/2-x(1/3-x(1/4-x(1/5-x(1/6)))))
c 
c        = 1 - x/2 + x^2/3 - x^3/4 + x^4/5 - x^5/6 
      implicit none
      real*4 x,log1pxox
      integer j,nterm
      if(x.le.-1.0) then
         stop 'Error: calling log1pxox with x < -1'
      elseif(abs(x).gt.0.1) then
         log1pxox=log(1+x)/x
      else
c         nterm=nint(2+110*abs(x))  ! double precision
         nterm=nint(2+50*abs(x))   ! single precision
         log1pxox=1.0/nterm
         do j=nterm-1,1,-1
            log1pxox=1.0/j-x*log1pxox
         end do
      endif
      return
      end
