      function log1pxox(x)
c  Computes log(1+x)/x without risk of underflow or zero-division. For x < 0.1
c  LOG1PXOX is approximated by the series 1-x/2+x**2/3-x**3/4+x**4/5-x**5/6
      real*4 x,log1pxox
      integer j,nterm
      if(abs(x).gt.0.1) then
         log1pxox=log(1+x)/x
      else
c         nterm=13  ! double precision
         nterm=6   ! single precision
         log1pxox=1.d0/nterm
         do j=nterm-1,1,-1
            log1pxox=1.d0/j-x*log1pxox
         end do
      endif
      return
      end
